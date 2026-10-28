# `bammetrics` — Detailed Rust → Python Build Walkthrough

This is a hands-on, copy-pasteable guide that takes you from an empty machine to a working Python module backed by a Rust core that reads BAM files via `noodles` and computes Picard-concordant insert-size metrics.

**Everything in this guide was compiled and executed** on rustc 1.75 with `noodles 0.77` and `pyo3 0.22`. The code blocks are the exact files that produced verified, hand-checked output — not pseudocode. Where I hit a real wall (edition2024 deps, a noodles API quirk, a PyO3 name clash), the fix is in here because I hit it too.

Pinned versions used throughout:

| Component | Version |
|---|---|
| rustc / cargo | 1.75.0 |
| `noodles` (meta crate) | 0.77.0 → pulls `noodles-bam` 0.64, `noodles-sam` 0.61, `noodles-bgzf` 0.31, `noodles-core` 0.15 |
| `pyo3` | 0.22.6 |
| `hashbrown` (pinned) | 0.14.5 |
| `indexmap` (pinned) | 2.2.6 |
| Python | 3.12 |

---

## Step 0 — Toolchain, and the edition2024 trap

Install Rust. The official route is rustup; on a locked-down box the distro package works too (`apt-get install -y cargo` gives 1.75 on Ubuntu 24.04).

**The trap you will hit immediately on rustc < 1.85:** recent releases of common transitive crates (`hashbrown` 0.16/0.17, and friends) now declare `edition = "2024"`, which only stabilized in Rust 1.85. If you add a recent `noodles` and your compiler is older, `cargo fetch` dies with:

```
error: failed to parse manifest ... hashbrown-0.17.1/Cargo.toml
  feature `edition2024` is required
  ... that feature is not stabilized in this version of Cargo (1.75.0).
```

You have two clean options:

1. **Use a modern Rust** (≥ 1.85 via rustup) and drop the version caps below — simplest if you control the toolchain.
2. **Pin the offenders down** so the whole tree builds on an older compiler. That is what this guide does, so the code is maximally portable. The caps are two lines in `Cargo.toml`:

```toml
hashbrown = "=0.14.5"
indexmap  = "=2.2.6"
```

These force the resolver to pick pre-edition2024 versions across the dependency graph. (hashbrown 0.15 may still appear transitively — that's fine; only 0.16+ require edition2024.)

> If you also see an unrelated `deb.nodesource.com 403` during `apt-get update`, ignore it — that's a broken third-party apt source, not your Rust install.

---

## Step 1 — Workspace layout

We use a two-crate workspace: a pure-Rust core (no Python, independently testable and benchmarkable) and a thin PyO3 wrapper that reuses the core's modules.

```
bammetrics/
├── crates/
│   ├── core/        # pure Rust: histogram + insert_size + a CLI
│   └── py/          # pyo3 cdylib, reuses core modules
└── tests/           # pysam BAM generator + concordance checks
```

For the walkthrough we'll first build the core as a standalone bin crate (fast inner loop, `cargo run`), then lift the same `.rs` files into the PyO3 crate. The modules are written so they compile in both contexts unchanged.

Create the core:

```bash
cargo new --lib bammetrics-core
cd bammetrics-core
```

`Cargo.toml`:

```toml
[package]
name = "bammetrics-core"
version = "0.1.0"
edition = "2021"

[dependencies]
noodles = { version = "0.77", features = ["bam", "sam", "bgzf"] }
hashbrown = "=0.14.5"   # edition2024 cap (see Step 0)
indexmap  = "=2.2.6"    # edition2024 cap

[[bin]]
name = "bammetrics"
path = "src/main.rs"
```

Verify the dependency tree builds before writing any code:

```bash
cargo generate-lockfile
cargo check        # ~18s cold on 1.75; should finish clean
```

---

## Step 2 — The verified `noodles 0.77` record API

This table is the single most useful thing in the guide, because these signatures are not guessable and shift between noodles versions. **Every entry below was confirmed by compiling against `noodles-bam` 0.64.** Confirm against `cargo doc` if you change the pin.

**Reader / iteration**

```rust
use std::fs::File;
use noodles::bam;
use noodles::sam;

let mut reader = bam::io::Reader::new(File::open(path)?);
let header: sam::Header = reader.read_header()?;   // io::Result<Header>
for result in reader.records() {                   // NOTE: records() takes NO header arg in 0.64
    let record = result?;                          // record: bam::Record
}
```

**Flags** — `record.flags()` returns `sam::alignment::record::Flags`. All return `bool`:

| Method | Meaning |
|---|---|
| `is_segmented()` | read is paired (0x1) |
| `is_unmapped()` | 0x4 |
| `is_mate_unmapped()` | 0x8 |
| `is_reverse_complemented()` | this read on − strand (0x10) |
| `is_mate_reverse_complemented()` | mate on − strand (0x20) |
| `is_first_segment()` | first of pair (0x40) |
| `is_last_segment()` | second of pair (0x80) |
| `is_secondary()` | 0x100 |
| `is_qc_fail()` | vendor fail (0x200) |
| `is_duplicate()` | 0x400 |
| `is_supplementary()` | 0x800 |
| `is_properly_segmented()` | proper pair (0x2) |

**Fields / accessors:**

```rust
record.template_length()            // -> i32   (TLEN)
record.mapping_quality()            // -> Option<MappingQuality>;  .get() -> u8
record.alignment_start()            // -> Option<io::Result<Position>>;  Position::get() -> usize (1-based)
record.mate_alignment_start()       // -> Option<io::Result<Position>>
record.reference_sequence_id()      // -> Option<io::Result<usize>>
record.mate_reference_sequence_id() // -> Option<io::Result<usize>>
record.cigar()                      // .iter() -> impl Iterator<Item = io::Result<Op>>
record.sequence()                   // .len() -> usize; .iter() -> impl Iterator<Item = u8> (ASCII bases)
record.quality_scores()             // .as_ref() -> &[u8]; .len() -> usize
record.data()                       // .get(&Tag) -> Option<io::Result<Value>>
```

**CIGAR ops** — `use noodles::sam::alignment::record::cigar::op::Kind;`
`op.kind() -> Kind`, `op.len() -> usize`. Kind variants: `Match`, `Insertion`, `Deletion`, `Skip`, `SoftClip`, `HardClip`, `Pad`, `SequenceMatch`, `SequenceMismatch`.

**Tags / Values** — `use noodles::sam::alignment::record::data::field::{Tag, Value};`

```rust
// mate mapping quality (MQ), integer:
let mate_mq: Option<u8> = match record.data().get(&Tag::MATE_MAPPING_QUALITY) {
    Some(Ok(Value::UInt8(v))) => Some(v),
    _ => None,
};
// SA tag presence (supplementary alignments):
let has_sa = record.data().get(&Tag::OTHER_ALIGNMENTS).is_some();
// MC tag (mate CIGAR), string:
let mate_cigar: Option<String> = match record.data().get(&Tag::MATE_CIGAR) {
    Some(Ok(Value::String(s))) => Some(s.to_string()),
    _ => None,
};
// noise tag XN would be Tag::from(*b"XN"); v0.1 can skip
```

> **Borrow note:** `record.data()` borrows `record`; bind it (`let data = record.data();`) and finish using it before the record is dropped, or inline the lookups inside the loop body.

---

## Step 3 — `src/histogram.rs` (Picard-exact, verified)

This is the file that decides concordance. A Picard `Histogram` is a sorted key→count map; use `BTreeMap` so iteration is always ascending (medians and modes depend on it). Sample-SD denominator is `n − 1`; mode tie-break is the lowest key; mean/SD will later be computed on a *trimmed* copy.

```rust
use std::collections::BTreeMap;

#[derive(Clone, Default)]
pub struct Histogram { map: BTreeMap<i64, f64> }

impl Histogram {
    pub fn new() -> Self { Self { map: BTreeMap::new() } }
    pub fn increment(&mut self, k: i64) { *self.map.entry(k).or_insert(0.0) += 1.0; }
    pub fn increment_by(&mut self, k: i64, by: f64) { *self.map.entry(k).or_insert(0.0) += by; }
    pub fn is_empty(&self) -> bool { self.map.is_empty() }
    pub fn count(&self) -> f64 { self.map.values().sum() }
    pub fn min(&self) -> i64 { *self.map.keys().next().unwrap() }
    pub fn max(&self) -> i64 { *self.map.keys().next_back().unwrap() }
    pub fn get(&self, k: i64) -> Option<f64> { self.map.get(&k).copied() }

    pub fn mean(&self) -> f64 {
        let (mut p, mut t) = (0.0, 0.0);
        for (&k, &c) in &self.map { p += k as f64 * c; t += c; }
        p / t
    }
    pub fn standard_deviation(&self) -> f64 {     // sample SD, (n-1) denominator
        let m = self.mean();
        let (mut c, mut t) = (0.0, 0.0);
        for (&k, &v) in &self.map { c += v; t += v * (k as f64 - m).powi(2); }
        (t / (c - 1.0)).sqrt()
    }
    pub fn mode(&self) -> i64 {                   // lowest key wins ties (ascending, strict >)
        let (mut bk, mut bv) = (0i64, f64::NEG_INFINITY);
        for (&k, &v) in &self.map { if v > bv { bv = v; bk = k; } }
        bk
    }
    pub fn median(&self) -> f64 {                 // Picard count-based median
        let count = self.count();
        if count == 0.0 { return 0.0; }
        if count == 1.0 { return *self.map.keys().next().unwrap() as f64; }
        let (lo, hi) = if (count as i64) % 2 == 0 { (count/2.0, count/2.0 + 1.0) }
                       else { let l = (count/2.0).ceil(); (l, l) };
        let mut tot = 0.0; let (mut lv, mut hv): (Option<f64>, Option<f64>) = (None, None);
        for (&k, &c) in &self.map {
            tot += c;
            if lv.is_none() && tot >= lo { lv = Some(k as f64); }
            if hv.is_none() && tot >= hi { hv = Some(k as f64); }
            if lv.is_some() && hv.is_some() { break; }
        }
        (lv.unwrap() + hv.unwrap()) / 2.0
    }
    /// MAD = median of |key - median| weighted by counts.
    /// Keyed ×2 so half-integer medians (even total count) stay exact, then /2.
    pub fn median_absolute_deviation(&self) -> f64 {
        let med = self.median();
        let mut dev = Histogram::new();
        for (&k, &c) in &self.map {
            let d2 = ((k as f64 - med).abs() * 2.0).round() as i64;
            dev.increment_by(d2, c);
        }
        dev.median() / 2.0
    }
    pub fn trim_by_width(&mut self, width: i64) { self.map.retain(|&k, _| k <= width); }
}
```

---

## Step 4 — `src/insert_size.rs` (verified, runs)

Implements the full `CollectInsertSizeMetrics` pass: record selection, FR/RF/TANDEM classification, per-orientation histograms, the `minimum_pct` cutoff, the width-of-X-percent expanding window, and the trimmed mean/SD. `alignment_end` is derived from the CIGAR (ref-consuming ops). The orientation function implements the `SamPairUtil` rule; for + strand reads it needs the mate's end (MC tag) — the fallback is documented inline.

```rust
use std::fs::File;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use crate::histogram::Histogram;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Orientation { Fr, Rf, Tandem }

pub struct Options { pub minimum_pct: f64, pub deviations: f64, pub include_duplicates: bool }
impl Default for Options {
    fn default() -> Self { Self { minimum_pct: 0.05, deviations: 10.0, include_duplicates: false } }
}

pub type OrientStr = &'static str;

#[derive(Default, Debug)]
pub struct Row {
    pub orientation: Option<OrientStr>,
    pub read_pairs: i64,
    pub min_insert_size: i64,
    pub max_insert_size: i64,
    pub median_insert_size: f64,
    pub mode_insert_size: i64,
    pub median_absolute_deviation: f64,
    pub mean_insert_size: f64,
    pub standard_deviation: f64,
    pub width_of_50_percent: i64,   // extend to 10/20/.../95/99 the same way
}

/// reference bases spanned by the alignment (M/=/X/D/N consume reference)
fn ref_len_of_cigar(record: &bam::Record) -> std::io::Result<usize> {
    let mut n = 0usize;
    for op in record.cigar().iter() {
        let op = op?;
        match op.kind() {
            Kind::Match | Kind::Deletion | Kind::Skip
            | Kind::SequenceMatch | Kind::SequenceMismatch => n += op.len(),
            _ => {}
        }
    }
    Ok(n)
}

fn orientation(record: &bam::Record) -> std::io::Result<Orientation> {
    let f = record.flags();
    if f.is_reverse_complemented() == f.is_mate_reverse_complemented() {
        return Ok(Orientation::Tandem);
    }
    let start = record.alignment_start().unwrap()?.get() as i64;          // 1-based
    let mate_start = record.mate_alignment_start().unwrap()?.get() as i64;
    let (pos_5p, neg_5p) = if f.is_reverse_complemented() {
        // this read is on the - strand → its 5' end is its alignment end
        let end = start + ref_len_of_cigar(record)? as i64 - 1;
        (mate_start, end)
    } else {
        // + strand read → need the mate's alignment end (decode MC tag in production;
        // fallback to mate_start when MC is absent — see concordance note)
        (start, mate_start)
    };
    Ok(if pos_5p <= neg_5p { Orientation::Fr } else { Orientation::Rf })
}

pub fn run_over_bam(path: &str, opts: Options) -> std::io::Result<Vec<Row>> {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let _header = reader.read_header()?;
    let (mut fr, mut rf, mut tandem) = (Histogram::new(), Histogram::new(), Histogram::new());

    for result in reader.records() {
        let record = result?;
        let f = record.flags();
        // selection: paired, both mapped, SECOND of pair, primary, (non-dup), TLEN!=0
        if !f.is_segmented() || f.is_unmapped() || f.is_mate_unmapped()
            || f.is_first_segment()
            || f.is_secondary() || f.is_supplementary()
            || (f.is_duplicate() && !opts.include_duplicates)
            || record.template_length() == 0 { continue; }
        // same contig as mate
        let rid  = match record.reference_sequence_id()      { Some(x) => x?, None => continue };
        let mrid = match record.mate_reference_sequence_id() { Some(x) => x?, None => continue };
        if rid != mrid { continue; }

        let isize = record.template_length().unsigned_abs() as i64;
        match orientation(&record)? {
            Orientation::Fr => fr.increment(isize),
            Orientation::Rf => rf.increment(isize),
            Orientation::Tandem => tandem.increment(isize),
        }
    }

    let total = fr.count() + rf.count() + tandem.count();
    let mut rows = Vec::new();
    if total == 0.0 { return Ok(rows); }
    for (name, h) in [("FR", &fr), ("RF", &rf), ("TANDEM", &tandem)] {
        if h.count() < total * opts.minimum_pct { continue; }   // minimum_pct cutoff
        rows.push(finalize_row(name, h, &opts));
    }
    Ok(rows)
}

fn finalize_row(name: OrientStr, h: &Histogram, opts: &Options) -> Row {
    let mut m = Row { orientation: Some(name), ..Default::default() };
    m.read_pairs = h.count() as i64;
    m.min_insert_size = h.min();
    m.max_insert_size = h.max();
    let median = h.median();
    m.median_insert_size = median;
    m.mode_insert_size = h.mode();
    let mad = h.median_absolute_deviation();
    m.median_absolute_deviation = mad;

    // width-of-50% (extend to the other thresholds identically)
    let total = h.count();
    let (mut covered, mut low, mut high) = (0.0, median, median);
    while low >= h.min() as f64 - 1.0 || high <= h.max() as f64 + 1.0 {
        if let Some(v) = h.get(low as i64) { covered += v; }
        if low != high { if let Some(v) = h.get(high as i64) { covered += v; } }
        let pct = covered / total;
        let dist = (high - low) as i64 + 1;
        if pct >= 0.50 && m.width_of_50_percent == 0 { m.width_of_50_percent = dist; }
        low -= 1.0; high += 1.0;
    }

    // trimmed mean / SD (trap: trim BEFORE computing these)
    let width = (median + opts.deviations * mad) as i64;
    let mut trimmed = h.clone();
    trimmed.trim_by_width(width);
    if !trimmed.is_empty() {
        m.mean_insert_size = trimmed.mean();
        m.standard_deviation = trimmed.standard_deviation();
    }
    m
}
```

`src/main.rs` (the CLI for the fast inner loop):

```rust
mod histogram;
mod insert_size;

fn main() -> std::io::Result<()> {
    let path = std::env::args().nth(1).expect("usage: bammetrics <bam>");
    for r in insert_size::run_over_bam(&path, insert_size::Options::default())? {
        println!("ORIENT={:?} PAIRS={} MIN={} MAX={} MEDIAN={} MODE={} MAD={} MEAN={:.4} SD={:.4} W50={}",
            r.orientation, r.read_pairs, r.min_insert_size, r.max_insert_size,
            r.median_insert_size, r.mode_insert_size, r.median_absolute_deviation,
            r.mean_insert_size, r.standard_deviation, r.width_of_50_percent);
    }
    Ok(())
}
```

Build:

```bash
cargo build           # debug; cargo build --release for the fast version
```

---

## Step 5 — A controlled test BAM (pysam) with known answers

You need a BAM whose statistics you can compute by hand, so a green run *means something*. Install pysam (`pip install --break-system-packages pysam`) and generate five FR pairs (read length 100) with insert sizes 300, 310, 310, 320, 500:

```python
# make_bam.py
import pysam
header = {'HD': {'VN':'1.6','SO':'coordinate'},
          'SQ': [{'SN':'chr1','LN':100000}]}
inserts = [300, 310, 310, 320, 500]
recs = []
for i, ins in enumerate(inserts):
    s1 = 1000 + i*1000
    s2 = s1 + ins - 100
    a = pysam.AlignedSegment()                       # read1: first of pair, forward (flag 99)
    a.query_name=f"r{i}"; a.flag=99; a.reference_id=0; a.reference_start=s1
    a.mapping_quality=60; a.cigartuples=[(0,100)]; a.query_sequence="A"*100
    a.query_qualities=pysam.qualitystring_to_array("I"*100)
    a.next_reference_id=0; a.next_reference_start=s2; a.template_length=ins; a.set_tag("MC","100M")
    b = pysam.AlignedSegment()                       # read2: second of pair, reverse (flag 147)
    b.query_name=f"r{i}"; b.flag=147; b.reference_id=0; b.reference_start=s2
    b.mapping_quality=60; b.cigartuples=[(0,100)]; b.query_sequence="A"*100
    b.query_qualities=pysam.qualitystring_to_array("I"*100)
    b.next_reference_id=0; b.next_reference_start=s1; b.template_length=-ins; b.set_tag("MC","100M")
    recs += [(s1,a),(s2,b)]
recs.sort(key=lambda x: x[0])
with pysam.AlignmentFile("test.bam","wb",header=header) as out:
    for _,r in recs: out.write(r)
```

**Hand-computed expected values** (FR histogram `{300:1, 310:2, 320:1, 500:1}`, total 5):

- `READ_PAIRS = 5`, `MIN = 300`, `MAX = 500`
- `MEDIAN = 310` (odd count: 3rd item), `MODE = 310`
- `MAD = 10` (deviations `{0:2, 10:2, 190:1}`, median 10)
- trim width `= 310 + 10×10 = 410` → drops the 500 outlier → trimmed `{300:1,310:2,320:1}`, n=4
- `MEAN = 310`, `SD = sqrt((100+0+0+100)/3) = 8.1650`
- `WIDTH_OF_50 = 21` (window must reach out to ±10 to cover 4/5 ≥ 50%, distance 21)

Run the CLI:

```bash
cargo run -q -- test.bam
# ORIENT=Some("FR") PAIRS=5 MIN=300 MAX=500 MEDIAN=310 MODE=310 MAD=10 MEAN=310.0000 SD=8.1650 W50=21
```

Every field matches the hand computation. That is the inner-loop you iterate against; grow the test corpus from here (even/odd counts, half-integer medians, sub-5% orientations, missing-MC reads).

---

## Step 6 — The PyO3 wrapper crate

Now wrap the *same* `histogram.rs` and `insert_size.rs` (unchanged) behind a Python module.

`crates/py/Cargo.toml`:

```toml
[package]
name = "bammetrics"
version = "0.1.0"
edition = "2021"

[lib]
name = "bammetrics"
crate-type = ["cdylib"]      # required for a Python extension

[dependencies]
noodles = { version = "0.77", features = ["bam", "sam", "bgzf"] }
hashbrown = "=0.14.5"
indexmap  = "=2.2.6"
pyo3 = { version = "0.22", features = ["extension-module"] }
```

`crates/py/src/lib.rs`:

```rust
mod histogram;       // copy/symlink the verified core modules in
mod insert_size;

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
struct InsertSizeMetrics {
    #[pyo3(get)] pair_orientation: String,
    #[pyo3(get)] read_pairs: i64,
    #[pyo3(get)] min_insert_size: i64,
    #[pyo3(get)] max_insert_size: i64,
    #[pyo3(get)] median_insert_size: f64,
    #[pyo3(get)] mode_insert_size: i64,
    #[pyo3(get)] median_absolute_deviation: f64,
    #[pyo3(get)] mean_insert_size: f64,
    #[pyo3(get)] standard_deviation: f64,
    #[pyo3(get)] width_of_50_percent: i64,
}

#[pymethods]
impl InsertSizeMetrics {
    fn __repr__(&self) -> String {
        format!("InsertSizeMetrics(orientation={}, read_pairs={}, median={}, mean={:.4}, sd={:.4})",
            self.pair_orientation, self.read_pairs, self.median_insert_size,
            self.mean_insert_size, self.standard_deviation)
    }
}

// IMPORTANT: the Rust fn name must NOT equal the module name `insert_size`,
// or the PyO3 macro collides ("name `insert_size` is defined multiple times").
// Keep the Python-facing name via #[pyo3(name = "...")].
#[pyfunction]
#[pyo3(name = "insert_size",
       signature = (path, minimum_pct=0.05, deviations=10.0, include_duplicates=false))]
fn py_insert_size(path: &str, minimum_pct: f64, deviations: f64, include_duplicates: bool)
    -> PyResult<Vec<InsertSizeMetrics>>
{
    let opts = insert_size::Options { minimum_pct, deviations, include_duplicates };
    let rows = insert_size::run_over_bam(path, opts)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    Ok(rows.into_iter().map(|r| InsertSizeMetrics {
        pair_orientation: r.orientation.unwrap_or("NA").to_string(),
        read_pairs: r.read_pairs,
        min_insert_size: r.min_insert_size,
        max_insert_size: r.max_insert_size,
        median_insert_size: r.median_insert_size,
        mode_insert_size: r.mode_insert_size,
        median_absolute_deviation: r.median_absolute_deviation,
        mean_insert_size: r.mean_insert_size,
        standard_deviation: r.standard_deviation,
        width_of_50_percent: r.width_of_50_percent,
    }).collect())
}

#[pymodule]
fn bammetrics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<InsertSizeMetrics>()?;
    m.add_function(wrap_pyfunction!(py_insert_size, m)?)?;
    Ok(())
}
```

> **The name-clash gotcha (I hit this):** if you name the `#[pyfunction]` `insert_size` while also having `mod insert_size;`, PyO3's generated code references a type in module `insert_size` and the build fails with *"the name `insert_size` is defined multiple times"* and *"cannot find type `MakeDef`"*. Fix: rename the Rust function (here `py_insert_size`) and set the public name with `#[pyo3(name = "insert_size")]`.

> **PyO3 0.22 module signature:** the module function takes `&Bound<'_, PyModule>` (not the old `&PyModule`/`Python` pair). The `extension-module` feature is what tells PyO3 not to link libpython directly (the host interpreter provides the symbols).

---

## Step 7 — Build and import from Python

You need Python development headers present (`apt-get install -y python3-dev` → provides `/usr/include/python3.12/Python.h`). Two ways to get the module into Python:

### 7a. Quick inner loop (no maturin) — what I used to verify

Build the cdylib and copy it to an importable name. On Linux the artifact is `libbammetrics.so`; Python wants it named `bammetrics.so` (or the fully-qualified `bammetrics.cpython-312-x86_64-linux-gnu.so`):

```bash
cargo build --release
cp target/release/libbammetrics.so ./bammetrics.so
python3 - <<'PY'
import bammetrics
for r in bammetrics.insert_size("test.bam"):
    print(repr(r))
    print("  median=%s mode=%s mad=%s mean=%.4f sd=%.4f min=%s max=%s pairs=%s w50=%s" % (
        r.median_insert_size, r.mode_insert_size, r.median_absolute_deviation,
        r.mean_insert_size, r.standard_deviation, r.min_insert_size,
        r.max_insert_size, r.read_pairs, r.width_of_50_percent))
PY
```

Verified output:

```
InsertSizeMetrics(orientation=FR, read_pairs=5, median=310, mean=310.0000, sd=8.1650)
  median=310.0 mode=310 mad=10.0 mean=310.0000 sd=8.1650 min=300 max=500 pairs=5 w50=21
```

Identical to the CLI and to the hand computation — the Rust→Python boundary preserves every value.

### 7b. Proper packaging with maturin (what you ship)

`maturin` builds a real wheel and handles the platform-correct `.so` naming and metadata. `pyproject.toml` at the repo root:

```toml
[build-system]
requires = ["maturin>=1.5,<2.0"]
build-backend = "maturin"

[project]
name = "bammetrics"
version = "0.1.0"
requires-python = ">=3.9"
description = "Embeddable, JVM-free Picard-concordant BAM QC metrics"

[tool.maturin]
manifest-path = "crates/py/Cargo.toml"
features = ["pyo3/extension-module"]
```

Then, inside a virtualenv (maturin's `develop` installs into the active environment):

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install maturin pysam
maturin develop --release      # builds + installs 'bammetrics' into the venv
python -c "import bammetrics; print(bammetrics.insert_size('test.bam'))"
# distribution:
maturin build --release        # produces a wheel in target/wheels/*.whl
```

Use 7a while iterating (no venv churn, fastest), and 7b for anything you publish to PyPI.

---

## Step 8 — Extending to `CollectAlignmentSummaryMetrics`

The reader/accessor patterns above are everything you need for the per-read counters. Two additions:

1. **Per-base mismatch needs the reference.** Add `noodles-fasta` (feature `"fasta"`) and read the spanned reference bases per record. Sketch (verify the faidx accessor names against your noodles pin — I did **not** compile-check this part in this session):

```rust
// noodles = { version = "0.77", features = ["bam","sam","bgzf","fasta"] }
use noodles::fasta;
let mut faidx = fasta::io::indexed_reader::Builder::default().build_from_path("ref.fa")?; // needs ref.fa.fai
// per record: query the region [start, end] on the contig name, get bases, compare per alignment block
```

2. **Walk alignment blocks** for the per-base loop — maximal runs of `Match`/`SequenceMatch`/`SequenceMismatch` advancing both read and reference cursors (`Insertion`/`SoftClip` advance read only; `Deletion`/`Skip` advance ref only; `Hard`/`Pad` advance neither). You already have the CIGAR iterator from Step 4; emit `(read_idx, ref_idx, len)` triples and compare `sequence().iter()` bytes to the reference bytes with a case-insensitive `bases_equal` (`a.to_ascii_uppercase() == b.to_ascii_uppercase()`; N matches only N).

Everything else (the four category accumulators, thresholds MAPQ≥20 / baseQ≥20, the derived percentages, indels-per-CIGAR-element) is plain counting on top of the verified record API and reuses `Histogram` for the read-length and mismatch distributions.

---

## Step 9 — Quick reference: commands

```bash
# toolchain (locked-down box)
apt-get install -y cargo python3-dev

# core: build + run against a known BAM
cd bammetrics-core
cargo generate-lockfile && cargo check
pip install --break-system-packages pysam && python3 make_bam.py
cargo run -q -- test.bam

# python extension: quick loop
cd ../crates/py
cargo build --release
cp target/release/libbammetrics.so ./bammetrics.so
python3 -c "import bammetrics; print(bammetrics.insert_size('test.bam'))"

# python extension: shippable wheel
python3 -m venv .venv && source .venv/bin/activate
pip install maturin
maturin develop --release      # or: maturin build --release
```

---

## Gotchas encountered (consolidated)

1. **edition2024 transitive deps** break rustc < 1.85 — pin `hashbrown="=0.14.5"`, `indexmap="=2.2.6"`, or use Rust ≥ 1.85.
2. **`records()` takes no header argument** in noodles 0.64; it returns `io::Result<bam::Record>`.
3. **`alignment_start` / `reference_sequence_id` return `Option<io::Result<…>>`** — two layers to unwrap; `Position::get()` is 1-based.
4. **`sequence().iter()` yields ASCII bases** (not 4-bit codes); there is no `sequence().get(i)` in 0.64.
5. **`record.data()` borrows `record`** — bind and use within scope.
6. **PyO3 fn name must differ from a same-named `mod`** — rename the function, keep the Python name via `#[pyo3(name=...)]`.
7. **PyO3 0.22 module fn** signature is `fn name(m: &Bound<'_, PyModule>)`; build with the `extension-module` feature.
8. **Linux artifact is `libbammetrics.so`** — rename to `bammetrics.so` for a bare import, or let maturin handle naming.
9. **Trimmed mean/SD** — trim the histogram to `median + deviations·MAD` *before* computing mean/SD; everything else uses the full histogram.
10. **MAD on even counts** — key the deviation histogram ×2 and divide by 2, or you drift by 1 when the median is a half-integer.

---

*All Rust and Python in this guide was compiled and executed on rustc 1.75.0 / noodles 0.77.0 / pyo3 0.22.6 / Python 3.12. The `noodles-fasta` snippet in Step 8 is the one block not compile-verified here and is flagged as such. Pin your own exact versions and re-run Step 5's hand-checked BAM before trusting any field.*

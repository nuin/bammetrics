# Synthetic corpus, Picard 3.5.0 references and the concordance harness.
# Needs samtools, bwa, picard (3.5.0) and uv on PATH.

UV_RUN := uv run --no-project
CORPUS := tests/corpus/out

.PHONY: all test lint release corpus reference-metrics concordance golden golden-update budgets a11y python-test wheel wheel-check clean-corpus bench-bam bench bench-wgs fuzz picard-fields
.PHONY: harness-test deny deflater-spike

all: test lint concordance

test:
	cargo test --workspace

lint:
	cargo fmt --all --check
	cargo clippy --workspace --all-targets -- -D warnings

release:
	cargo build --release

# Reference, intervals, t0-t3 BAMs. Rebuilt only when the generator changes.
corpus: $(CORPUS)/corpus.json

$(CORPUS)/corpus.json: tests/corpus/generate.py
	$(UV_RUN) --with pysam --with numpy python tests/corpus/generate.py $(CORPUS)

# Picard output for every case in tests/concordance/cases.py.
reference-metrics: corpus
	$(UV_RUN) python tests/concordance/run_picard.py

concordance: reference-metrics release golden harness-test
	$(UV_RUN) --with pytest --with pysam python -m pytest tests/concordance -q

# The harness's own tests (tests/harness): the Picard queue, compare_bam and
# compare_vcf on stand-in tools, the spec script. Needs samtools; starts no JVM.
harness-test:
	$(UV_RUN) --with pytest python -m pytest tests/harness -q

# Golden suite (crates/cli/tests/golden.rs): every documented command form on
# tests/data (and on corpus t1 with GOLDEN_CORPUS=1, since those records depend
# on one machine's bwa output), byte for byte
# against tests/golden. golden-update rewrites the records after an intended
# change; the diff is the review.
golden: release
	cargo test --release -p bammetrics-cli --test golden

golden-update: release
	UPDATE_GOLDEN=1 cargo test --release -p bammetrics-cli --test golden

# Size and speed budgets of docs/05 section 4.8 (crates/cli/tests/budgets.rs):
# the fixture pages and a 100-sample index timed (median of 31 renders) and
# sized, the binary growth between the two probe examples of crates/report
# and between the two of crates/formats, and the BAM reads the report adds to
# qc. Fails at 2x any target.
budgets:
	cargo build --release -p bammetrics-report --example size_base --example size_report
	cargo build --release -p bammetrics-formats --example size_formats_base --example size_formats
	cargo test --release -p bammetrics-cli --test budgets -- --ignored --nocapture

# Accessibility check of the report pages (benches/a11y.mjs): axe-core in
# Chromium, light and dark, 320, 360 and 1280 px, no sideways scroll at phone
# width with details closed and open, every control reached by Tab. Needs node
# and npm; installs Playwright and axe-core under target/a11y. Not run in CI.
A11Y := target/a11y
a11y: release
	npm install --no-audit --no-fund --prefix $(A11Y) playwright@1.64.1 axe-core@4.10.2
	cd $(A11Y) && npx playwright install chromium
	rm -rf $(A11Y)/pages $(A11Y)/run $(A11Y)/qc && mkdir -p $(A11Y)/run
	cargo run --release -q -p bammetrics-report --example fixture_pages -- $(A11Y)/pages
	cp -R tests/data/picard/ $(A11Y)/run/
	target/release/bammetrics report $(A11Y)/run --preset clinical-panel --date 2026-01-01 > /dev/null
	target/release/bammetrics qc tests/data/insert_size.bam -o $(A11Y)/qc/S1 --preset clinical-panel --panel tests/data/targets.interval_list --sex male --date "2026-01-01 00:00:00 UTC" > /dev/null
	node benches/a11y.mjs $(A11Y) tests/report/index.html $(A11Y)/pages/*.html $(A11Y)/run/*.html $(A11Y)/qc/*.html

# Differential fuzzing (tests/concordance/fuzz.py): N seeded random BAMs through
# Picard and bammetrics, failures shrunk into tests/concordance/regressions.
SEED ?= 1
N ?= 200
JOBS ?= 8
fuzz: release
	$(UV_RUN) python tests/concordance/fuzz.py --seed $(SEED) --n $(N) --jobs $(JOBS)

# The Picard metrics classes CompareMetrics reads (crates/core/picard_fields.tsv):
# every MetricBase subclass in picard.jar with its public fields, types and
# defaults, by reflection (tests/data/PicardFields.java). Rerun on a Picard upgrade.
PICARD_JAR ?= /opt/homebrew/Cellar/picard-tools/3.5.0/libexec/picard.jar
JAVA ?= /opt/homebrew/opt/openjdk/bin/java
picard-fields:
	$(JAVA) -cp $(PICARD_JAR) tests/data/PicardFields.java $(PICARD_JAR) > crates/core/picard_fields.tsv

# Python module tests: build into .venv with maturin, then run tests/python.
# Needs .venv with maturin, pytest and pysam, plus `make corpus` and `make release`.
python-test: corpus release
	.venv/bin/maturin develop --release
	.venv/bin/pytest tests/python -q

wheel:
	.venv/bin/maturin build --release

# The wheel as a user installs it: maturin build, a fresh venv with only the
# wheel and pytest, an import that must resolve inside that venv, then the
# tests that compare bammetrics.pyi with dir(bammetrics) both ways, and the
# stub the wheel ships against bammetrics.pyi, then mypy over the stub and the
# calls in tests/python/stub_types.py. Catches a stale wheel or stub.
WHEEL_CHECK := target/wheel-check
wheel-check:
	rm -rf $(WHEEL_CHECK)
	.venv/bin/maturin build --release -i python3 --out $(WHEEL_CHECK)/dist
	python3 -m venv $(WHEEL_CHECK)/venv
	$(WHEEL_CHECK)/venv/bin/pip install -q pytest mypy $(WHEEL_CHECK)/dist/*.whl
	$(WHEEL_CHECK)/venv/bin/mypy --cache-dir=$(WHEEL_CHECK)/mypy-cache bammetrics.pyi tests/python/stub_types.py
	cd $(WHEEL_CHECK) && venv/bin/python -c "import bammetrics, sys; f = bammetrics.__file__; \
		assert 'wheel-check/venv' in f, f; print('bammetrics', f)"
	cmp bammetrics.pyi $(WHEEL_CHECK)/venv/lib/python*/site-packages/bammetrics/__init__.pyi
	$(WHEEL_CHECK)/venv/bin/pytest -q -p no:cacheprovider tests/python/test_mirror.py \
		tests/python/test_api.py -k "stub or all_lists"

# 2M-pair benchmark BAM in target/bench (about 10 minutes, 320 MB; not committed).
bench-bam: target/bench/bench.json

target/bench/bench.json: benches/make_bam.py tests/corpus/generate.py
	$(UV_RUN) --with pysam --with numpy python benches/make_bam.py

# Picard vs bammetrics on the benchmark BAM; numbers go into docs/benchmarks.md.
bench: bench-bam release
	$(UV_RUN) python benches/bench.py

# T4 scale run of the coverage tools: 2 Mb reference, 200k pairs in target/bench_wgs
# (built on first use), Picard vs bammetrics; numbers go into docs/benchmarks.md.
bench-wgs: release
	$(UV_RUN) --with pysam --with numpy python benches/wgs_bench.py

clean-corpus:
	rm -rf $(CORPUS) target/concordance target/fakeR target/picard-queue

# Licences, banned and duplicated crates, registries and RustSec advisories of
# the dependency graph (deny.toml). Needs cargo-deny: cargo install cargo-deny --locked.
# Not run in CI.
deny:
	cargo deny check

# Deflater spike of docs/concordance.md ("BAM bytes"): Picard MarkDuplicates with
# each of its deflaters against bammetrics and three Rust-callable deflaters
# (benches/deflater_spike, a package outside the workspace). Needs picard.
SPIKE_BAMS ?= tests/data/markdup.bam
deflater-spike: release
	python3 benches/deflater_spike.py $(SPIKE_BAMS)

# Synthetic corpus, Picard 3.5.0 references and the concordance harness.
# Needs samtools, bwa, picard (3.5.0) and uv on PATH.

UV_RUN := uv run --no-project
CORPUS := tests/corpus/out

.PHONY: all test lint release corpus reference-metrics concordance golden golden-update python-test wheel wheel-check clean-corpus bench-bam bench bench-wgs fuzz

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

concordance: reference-metrics release golden
	$(UV_RUN) --with pytest --with pysam python -m pytest tests/concordance -q

# Golden suite (crates/cli/tests/golden.rs): every documented command form on
# tests/data (and on corpus t1 with GOLDEN_CORPUS=1, since those records depend
# on one machine's bwa output), byte for byte
# against tests/golden. golden-update rewrites the records after an intended
# change; the diff is the review.
golden: release
	cargo test --release -p bammetrics-cli --test golden

golden-update: release
	UPDATE_GOLDEN=1 cargo test --release -p bammetrics-cli --test golden

# Differential fuzzing (tests/concordance/fuzz.py): N seeded random BAMs through
# Picard and bammetrics, failures shrunk into tests/concordance/regressions.
SEED ?= 1
N ?= 200
JOBS ?= 8
fuzz: release
	$(UV_RUN) python tests/concordance/fuzz.py --seed $(SEED) --n $(N) --jobs $(JOBS)

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
	rm -rf $(CORPUS) target/concordance target/fakeR

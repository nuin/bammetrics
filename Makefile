# Synthetic corpus, Picard 3.5.0 references and the concordance harness.
# Needs samtools, bwa, picard (3.5.0) and uv on PATH.

UV_RUN := uv run --no-project
CORPUS := tests/corpus/out

.PHONY: all test lint release corpus reference-metrics concordance python-test wheel clean-corpus bench-bam bench

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

concordance: reference-metrics release
	$(UV_RUN) --with pytest python -m pytest tests/concordance -q

# Python module tests: build into .venv with maturin, then run tests/python.
# Needs .venv with maturin, pytest and pysam, plus `make corpus` and `make release`.
python-test: corpus release
	.venv/bin/maturin develop --release
	.venv/bin/pytest tests/python -q

wheel:
	.venv/bin/maturin build --release

# 2M-pair benchmark BAM in target/bench (about 10 minutes, 320 MB; not committed).
bench-bam: target/bench/bench.json

target/bench/bench.json: benches/make_bam.py tests/corpus/generate.py
	$(UV_RUN) --with pysam --with numpy python benches/make_bam.py

# Picard vs bammetrics on the benchmark BAM; numbers go into docs/benchmarks.md.
bench: bench-bam release
	$(UV_RUN) python benches/bench.py

clean-corpus:
	rm -rf $(CORPUS) target/concordance target/fakeR

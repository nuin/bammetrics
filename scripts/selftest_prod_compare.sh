#!/usr/bin/env bash
# End-to-end check of prod_compare.py on a fake pipeline run built from
# synthetic data. Picard stands in for bammetrics, so every metric must come
# back "match"; any other status is a bug in the comparison script.
#
# Needs picard, bwa, samtools and uv on PATH. Writes under target/selftest/.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
work="$repo/target/selftest"
rm -rf "$work"
run="$work/runs/260101_FAKE/Fake_2026_NGS_1"
mkdir -p "$run/BAM" "$work/ref" "$work/bin"

uv run --no-project --with pysam python - "$work" <<'EOF'
import random, sys
work = sys.argv[1]
rng = random.Random(7)
seq = "".join(rng.choice("ACGT") for _ in range(60000))
with open(f"{work}/ref/ref.fa", "w") as f:
    f.write(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
targets = [(2000 + i * 5000, 2000 + i * 5000 + 300) for i in range(10)]
with open(f"{work}/ref/panel.bed", "w") as f:
    for s, e in targets[:6]:
        f.write(f"chr1\t{s}\t{e}\tGENE{s}\n")
with open(f"{work}/ref/bait.bed", "w") as f:
    for s, e in targets:
        f.write(f"chr1\t{s - 50}\t{e + 50}\tBAIT{s}\n")
comp = str.maketrans("ACGT", "TGCA")
for sample in ("S1", "S2"):
    with open(f"{work}/{sample}_1.fq", "w") as r1, open(f"{work}/{sample}_2.fq", "w") as r2:
        for n in range(3000):
            s, e = rng.choice(targets) if rng.random() < 0.8 else (0, 59000)
            start = rng.randint(max(0, s - 200), min(e, 59400))
            frag = seq[start:start + rng.randint(250, 400)]
            a = list(frag[:100]); b = list(frag[-100:].translate(comp)[::-1])
            for read in (a, b):
                for i in range(len(read)):
                    if rng.random() < 0.005:
                        read[i] = rng.choice("ACGT")
            q = "".join(chr(33 + rng.choice((2, 15, 30, 37, 37, 37))) for _ in range(100))
            r1.write(f"@{sample}_{n}\n{''.join(a)}\n+\n{q}\n")
            r2.write(f"@{sample}_{n}\n{''.join(b)}\n+\n{q}\n")
EOF

ref="$work/ref/ref.fa"
samtools faidx "$ref"
picard CreateSequenceDictionary R="$ref" O="$work/ref/ref.dict" QUIET=true 2>/dev/null
bwa index "$ref" 2>/dev/null
for bed in panel bait; do
  picard BedToIntervalList I="$work/ref/$bed.bed" O="$work/ref/$bed.picard.bed" \
    SD="$work/ref/ref.dict" QUIET=true 2>/dev/null
done

for s in S1 S2; do
  mkdir -p "$run/BAM/$s/BAM" "$run/BAM/$s/Metrics"
  bwa mem -R "@RG\tID:$s\tSM:$s\tLB:$s\tPL:ILLUMINA" "$ref" "$work/${s}_1.fq" "$work/${s}_2.fq" 2>/dev/null |
    samtools sort -o "$run/BAM/$s/BAM/$s.bam" - 2>/dev/null
  samtools index "$run/BAM/$s/BAM/$s.bam"
done

cat > "$run/Fake.yaml" <<EOF
Reference: $ref
VCF: none
BAIT: $work/ref/bait.picard.bed
BED:
  S1: $work/ref/panel.bed
  S2: $work/ref/panel.bed
EOF

# Picard as the "existing" pipeline outputs under Metrics/, then as the
# stand-in bammetrics binary.
printf '#!/bin/sh\nexec picard "$@"\n' > "$work/bin/fake-bammetrics"
chmod +x "$work/bin/fake-bammetrics"
uv run --no-project --with pyyaml python "$repo/scripts/prod_compare.py" --datadir "$run" --bammetrics picard \
  --out "$work/seed" --metrics yield,align,hs,hs_panel,coverage,coverage_panel >/dev/null
for s in S1 S2; do cp "$work"/seed/$s/bammetrics/* "$run/BAM/$s/Metrics/"; done

uv run --no-project --with pyyaml python "$repo/scripts/prod_compare.py" \
  --datadir "$run" --bammetrics "$work/bin/fake-bammetrics" --picard picard \
  --out "$work/compare" --rerun-picard | tee "$work/summary.txt"

if grep -E '^S[12] ' "$work/summary.txt" | grep -vq ' match '; then
  echo "SELFTEST FAILED: a metric did not match" >&2
  exit 1
fi
echo "SELFTEST OK"

#!/usr/bin/env bash
# Wall time of one musl and one glibc bammetrics binary on the same commit.
#
#   benches/musl_speed.sh ARCH MUSL_BIN GLIBC_BIN
#
# ARCH is arm64 or amd64 (amd64 runs emulated under Docker Desktop on Apple silicon).
# Six commands on target/bench_wgs (make bench-wgs): the four commands of
# benches/wgs_bench.py, MarkDuplicates and CollectMultipleMetrics. Each runs 5
# times in one container per libc; wall time from /usr/bin/time, median of 5.
# Prints a Markdown table with musl/glibc ratios and their geometric mean.
# The musl binary runs in RUNNER_IMAGE (alpine), the glibc binary in
# RUST_GLIBC_IMAGE (bookworm), both from scripts/release.env.
# Environment: REPS (default 5), BENCH_DIR (default target/bench_wgs).
set -euo pipefail

[ "$#" = 3 ] || { echo "usage: $0 ARCH MUSL_BIN GLIBC_BIN" >&2; exit 2; }
arch="$1"
case "$arch" in arm64 | amd64) ;; *) echo "$0: ARCH is arm64 or amd64" >&2; exit 2 ;; esac
repo="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=../scripts/release.env
. "$repo/scripts/release.env"
musl_bin="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
glibc_bin="$(cd "$(dirname "$3")" && pwd)/$(basename "$3")"
bench="$(cd "${BENCH_DIR:-$repo/target/bench_wgs}" && pwd)"
reps="${REPS:-5}"
[ -f "$bench/bench.bam" ] || { echo "$0: $bench/bench.bam missing, run make bench-wgs" >&2; exit 2; }

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

inner='
set -eu
bam=/bench/bench.bam ref=/bench/reference.fa
silent="VALIDATION_STRINGENCY=SILENT VERBOSITY=ERROR QUIET=true"
measure() {
  label="$1"; shift
  : > /out/t.$$
  i=0
  while [ "$i" -lt "$REPS" ]; do
    /usr/bin/time -f %e -o /out/one "$@" > /dev/null 2> /out/err || { cat /out/err >&2; exit 1; }
    cat /out/one >> /out/t.$$
    i=$((i + 1))
  done
  m=$(sort -n /out/t.$$ | sed -n "$(( (REPS + 1) / 2 ))p")
  echo "$label $m"
  rm -f /out/t.$$ /out/*.txt /out/mm.* /out/md.bam
}
measure CollectWgsMetrics /bm CollectWgsMetrics I=$bam O=/out/a.txt R=$ref
measure "CollectWgsMetrics USE_FAST_ALGORITHM=true" /bm CollectWgsMetrics I=$bam O=/out/a.txt R=$ref USE_FAST_ALGORITHM=true
measure CollectRawWgsMetrics /bm CollectRawWgsMetrics I=$bam O=/out/a.txt R=$ref
measure CollectWgsMetricsWithNonZeroCoverage /bm CollectWgsMetricsWithNonZeroCoverage I=$bam O=/out/a.txt R=$ref
measure MarkDuplicates /bm MarkDuplicates I=$bam O=/out/md.bam M=/out/a.txt
measure CollectMultipleMetrics /bm CollectMultipleMetrics I=$bam O=/out/mm R=$ref
'

run_side() { # side image binary install-command
  local side="$1" image="$2" bin="$3" install="$4" out="$work/$1"
  mkdir -p "$out"
  docker run --rm --platform "linux/$arch" \
    -v "$bin:/bm:ro" -v "$bench:/bench:ro" -v "$out:/out" \
    -e REPS="$reps" "$image" sh -euc "$install; $inner" > "$work/$side.txt"
}
run_side musl "$RUNNER_IMAGE" "$musl_bin" "apk add --no-cache time >/dev/null"
run_side glibc "$RUST_GLIBC_IMAGE" "$glibc_bin" "apt-get update -qq >/dev/null && apt-get install -y -qq time >/dev/null"

python3 - "$work/glibc.txt" "$work/musl.txt" "$arch" "$reps" <<'PY'
import math, sys
def read(p):
    rows = []
    for line in open(p):
        label, _, t = line.rstrip("\n").rpartition(" ")
        rows.append((label, float(t)))
    return rows
g, m = read(sys.argv[1]), read(sys.argv[2])
assert [a for a, _ in g] == [a for a, _ in m], "command lists differ"
print(f"{sys.argv[3]}, median of {sys.argv[4]} wall times\n")
print("| command | glibc s | musl s | musl/glibc |")
print("|---|---:|---:|---:|")
ratios = []
for (label, gt), (_, mt) in zip(g, m):
    r = mt / gt
    ratios.append(r)
    print(f"| {label} | {gt:.2f} | {mt:.2f} | {r:.2f} |")
mean = math.exp(sum(map(math.log, ratios)) / len(ratios))
print(f"| geometric mean | | | {mean:.2f} |")
PY

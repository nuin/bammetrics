#!/usr/bin/env bash
# Build a static x86_64 Linux binary of the bammetrics CLI for prod3/prod4.
#
# The hosts run an old glibc and have no internet, so the binary is linked
# against musl inside an amd64 Alpine container: it has no shared-library
# dependencies and runs on any x86_64 Linux kernel. On an Apple Silicon Mac
# Docker runs the amd64 container under emulation, so the first build takes
# several minutes; the cargo registry and target dir are cached in volumes.
#
# Output: dist/bammetrics-x86_64-linux-musl and its .sha256.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
image="rust:1-alpine"
target="x86_64-unknown-linux-musl"

manifest="crates/core/Cargo.toml"
if grep -q '^\[workspace\]' "$repo/Cargo.toml" 2>/dev/null; then
  manifest="Cargo.toml"
fi

docker run --rm --platform linux/amd64 \
  -v "$repo":/src -w /src \
  -v bammetrics-cargo-registry:/usr/local/cargo/registry \
  -v bammetrics-linux-target:/src/target-linux \
  -e CARGO_TARGET_DIR=/src/target-linux \
  "$image" sh -euc "
    apk add --no-cache musl-dev >/dev/null
    cargo build --release --locked --manifest-path $manifest --bin bammetrics --target $target
    mkdir -p dist
    cp target-linux/$target/release/bammetrics dist/bammetrics-x86_64-linux-musl
  "

cd "$repo/dist"
shasum -a 256 bammetrics-x86_64-linux-musl > bammetrics-x86_64-linux-musl.sha256
file bammetrics-x86_64-linux-musl
cat bammetrics-x86_64-linux-musl.sha256

#!/usr/bin/env bash
# Pack the Linux binary, its checksum, the comparison script and its doc into
# dist/bammetrics-prodtest.tar.gz for copying to prod3/prod4.
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo"
test -f dist/bammetrics-x86_64-linux-musl || { echo "run scripts/build-linux.sh first" >&2; exit 1; }
stage="$(mktemp -d)"
cp dist/bammetrics-x86_64-linux-musl dist/bammetrics-x86_64-linux-musl.sha256 \
   scripts/prod_compare.py docs/prod-validation.md "$stage/"
tar -czf dist/bammetrics-prodtest.tar.gz -C "$stage" .
rm -rf "$stage"
tar -tzf dist/bammetrics-prodtest.tar.gz
ls -la dist/bammetrics-prodtest.tar.gz

#!/usr/bin/env bash
# Build the BAMmetrics release artifacts from one commit: two static musl
# Linux binaries (Docker, pinned images), the macOS universal binary (rustup
# toolchain pinned in release.env), one tarball each, SHA256SUMS, BUILDINFO.
#
#   scripts/release.sh [--snapshot | --tag] [--dry-run] [--publish] [--clean]
#
# --snapshot (default)  dist/<version>-<shortsha>/, any commit
# --tag                 dist/v<version>/, HEAD must carry the tag v<version>
# --dry-run             print one line per command, run none
# --publish             gh release create as a draft, then published (needs --tag)
# --clean               drop the build volumes and target/release-mac first
#
# Exit 2: refused before any build. Exit 1: a build or a self-check failed.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo"
blocked_from_env="${BLOCKED_TARGETS-}"
# shellcheck source=release.env
. "$repo/scripts/release.env"
BLOCKED_TARGETS="${BLOCKED_TARGETS-} $blocked_from_env"

mode=snapshot dry=0 publish=0 clean=0
die() { echo "release.sh: $*" >&2; exit 2; }
for arg in "$@"; do
  case "$arg" in
    --snapshot) mode=snapshot ;;
    --tag) mode=tag ;;
    --dry-run) dry=1 ;;
    --publish) publish=1 ;;
    --clean) clean=1 ;;
    *) die "unknown option $arg (--snapshot, --tag, --dry-run, --publish, --clean)" ;;
  esac
done
[ "$publish" = 0 ] || [ "$mode" = tag ] || die "--publish needs --tag"

# --- refusals, before any build -------------------------------------------
dirty="$( { git diff --name-only HEAD; git diff --cached --name-only; } | sort -u)"
[ -z "$dirty" ] || die "tracked files differ from HEAD: $(echo "$dirty" | tr '\n' ' ')"

untracked="$(git ls-files --others --exclude-standard -- .cargo crates vendor)"
[ -z "$untracked" ] || die "untracked files under .cargo, crates or vendor: $(echo "$untracked" | tr '\n' ' ')"

version="$(grep -m1 '^version *= *"' Cargo.toml | cut -d'"' -f2)"
[ -n "$version" ] || die "no version in Cargo.toml"
commit="$(git rev-parse --short HEAD)"
full_commit="$(git rev-parse HEAD)"
if [ "$mode" = tag ]; then
  tag_commit="$(git rev-parse -q --verify "refs/tags/v$version^{commit}" || true)"
  [ "$tag_commit" = "$full_commit" ] || die "tag v$version does not point at HEAD ($commit)"
  out="dist/v$version"
else
  out="dist/$version-$commit"
fi

notes_text=""
if [ "$publish" = 1 ]; then
  notes_text="$(awk -v h="## $version" '$0 == h {on=1; next} /^## / {on=0} on' CHANGELOG.md)"
  [ -n "$(printf '%s' "$notes_text" | tr -d '[:space:]')" ] \
    || die "CHANGELOG.md has no notes under '## $version'"
  # The tag on origin must be the commit being released, or gh would tag another one.
  if [ "$dry" = 1 ]; then
    echo "+ git ls-remote --tags origin refs/tags/v$version (must resolve to $full_commit)"
  else
    remote_tag="$(git ls-remote --tags origin "refs/tags/v$version" "refs/tags/v$version^{}" | awk '{print $1}' | tail -n 1)"
    [ "$remote_tag" = "$full_commit" ] \
      || die "tag v$version on origin is '${remote_tag:-missing}', not $full_commit"
  fi
fi

need() { command -v "$1" >/dev/null 2>&1 || die "missing tool: $1"; }
for tool in docker lipo shasum python3 file arch; do need "$tool"; done
[ "$publish" = 0 ] || need gh
rustup="$(command -v rustup || true)"
if [ -z "$rustup" ] && command -v brew >/dev/null 2>&1; then
  candidate="$(brew --prefix rustup 2>/dev/null)/bin/rustup"
  [ -x "$candidate" ] && rustup="$candidate"
fi
[ -n "$rustup" ] || die "missing tool: rustup"
installed="$("$rustup" target list --toolchain "$MAC_TOOLCHAIN" --installed 2>/dev/null || true)"
for t in aarch64-apple-darwin x86_64-apple-darwin; do
  echo "$installed" | grep -qx "$t" || die "rustup toolchain $MAC_TOOLCHAIN lacks target $t"
done

# --- helpers ----------------------------------------------------------------
SOURCE_DATE_EPOCH="$(git log -1 --format=%ct)"
export SOURCE_DATE_EPOCH
cargo_home="${CARGO_HOME:-$HOME/.cargo}"
mac_rustflags="--remap-path-prefix=$repo=/build --remap-path-prefix=$cargo_home=/cargo"
linux_rustflags="--remap-path-prefix=$repo=/build --remap-path-prefix=/usr/local/cargo=/cargo --remap-path-prefix=/target=/build/target"
git_common="$(cd "$repo" && cd "$(git rev-parse --git-common-dir)" && pwd)"
work="$out/.work"

run() {
  if [ "$dry" = 1 ]; then printf '+'; printf ' %s' "$@"; echo; else "$@"; fi
}
# Shell text for a container script: install NAME@VERSION unless that version is
# in the tools volume, then fail when it still is not.
pin_tool() {
  printf "grep -q '^\"%s %s (' /tools/.crates.toml 2>/dev/null || cargo install --locked %s@%s; grep -q '^\"%s %s (' /tools/.crates.toml || { echo 'tool %s %s is not installed' >&2; exit 1; }" \
    "$1" "$2" "$1" "$2" "$1" "$2" "$1" "$2"
}
blocked() { case " $BLOCKED_TARGETS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
fail() { echo "release.sh: self-check failed: $*" >&2; exit 1; }

targets=()
for t in x86_64-unknown-linux-musl aarch64-unknown-linux-musl universal-apple-darwin; do
  if blocked "$t"; then
    echo "release.sh: skipping $t (blocked in release.env)" >&2
  else
    targets+=("$t")
  fi
done
has() { local t; for t in "${targets[@]}"; do [ "$t" = "$1" ] && return 0; done; return 1; }
arch_of() { case "$1" in x86_64-*) echo amd64 ;; *) echo arm64 ;; esac; }
name_of() { echo "bammetrics-$version-$1"; }

# --- builds ----------------------------------------------------------------
if [ "$clean" = 1 ]; then
  for t in x86_64-unknown-linux-musl aarch64-unknown-linux-musl; do
    run docker volume rm -f "bammetrics-release-target-$(arch_of "$t")"
  done
  run rm -rf "$repo/target/release-mac" "$repo/$out"
fi
# Files left by an earlier run in this directory are not part of this run, and a
# run that fails must not leave old sums beside new files.
if [ "$dry" = 1 ]; then
  echo "+ rm -f $out/*.tar.gz $out/*.cdx.json $out/bammetrics.rb $out/SHA256SUMS $out/BUILDINFO"
else
  mkdir -p "$out"
  rm -f "$out"/*.tar.gz "$out"/*.cdx.json "$out"/bammetrics.rb "$out"/SHA256SUMS "$out"/BUILDINFO
fi
[ "$dry" = 1 ] || mkdir -p "$work"

linux_build() {
  local target="$1" arch bin_dir script
  arch="$(arch_of "$target")"
  bin_dir="$repo/$work/$(name_of "$target")"
  script="apk add --no-cache musl-dev git >/dev/null; git config --global --add safe.directory '*'"
  script+="; $(pin_tool cargo-auditable "$CARGO_AUDITABLE_VERSION")"
  script+="; apk list -I musl-dev | cut -d\" \" -f1 > $repo/$work/musl-dev.txt"
  script+="; export PATH=/tools/bin:\$PATH"
  script+="; cargo auditable build --release --locked -p bammetrics-cli --features pdf --target $target"
  script+="; rustc --version > $repo/$work/rustc-$target.txt"
  script+="; mkdir -p $bin_dir; cp /target/$target/release/bammetrics $bin_dir/bammetrics"
  script+="; chown -R $(id -u):$(id -g) $repo/$work"
  run docker run --rm --platform "linux/$arch" \
    -v "$repo:$repo" -w "$repo" \
    -v "$git_common:$git_common:ro" \
    -v bammetrics-cargo-registry:/usr/local/cargo/registry \
    -v bammetrics-cargo-tools:/tools \
    -v "bammetrics-release-target-$arch:/target" \
    -e CARGO_TARGET_DIR=/target -e CARGO_INSTALL_ROOT=/tools \
    -e SOURCE_DATE_EPOCH="$SOURCE_DATE_EPOCH" -e RUSTFLAGS="$linux_rustflags" \
    "$RUST_ALPINE_IMAGE" sh -euc "$script"
}

linux_check() {
  local target="$1" arch bin got
  arch="$(arch_of "$target")"
  bin="$repo/$work/$(name_of "$target")/bammetrics"
  if [ "$dry" = 1 ]; then
    echo "+ self-check $target: --version --json, file"
    return
  fi
  got="$(BM_LINUX_BIN="$bin" BM_LINUX_ARCH="$arch" "$repo/scripts/linux-bammetrics" --version --json)" \
    || fail "$target does not run"
  version_check "$target" "$got"
  file "$bin" | grep -Eq 'ELF 64-bit.*(statically linked|static-pie linked)' \
    || fail "$target is not a static ELF: $(file "$bin")"
  file "$bin" | grep -q "$([ "$arch" = amd64 ] && echo x86-64 || echo aarch64)" \
    || fail "$target has the wrong architecture: $(file "$bin")"
  echo "self-check $target: $got"
}

version_check() {
  python3 -c '
import json, sys
d = json.loads(sys.argv[3])
ok = d.get("version") == sys.argv[1] and d.get("commit") == sys.argv[2]
sys.exit(0 if ok else 1)' "$version" "$commit" "$2" \
    || fail "$1 reports $2, expected version $version commit $commit"
}

mac_build() {
  local target="$1" sysroot="<sysroot>"
  # Homebrew's cargo can come first on PATH even under `rustup run`; put the
  # pinned toolchain's bin in front so its rustc, with both Apple targets, builds.
  [ "$dry" = 1 ] || sysroot="$("$rustup" run "$MAC_TOOLCHAIN" rustc --print sysroot)"
  run env PATH="$sysroot/bin:$PATH" RUSTC="$sysroot/bin/rustc" CARGO_TARGET_DIR="$repo/target/release-mac" \
    MACOSX_DEPLOYMENT_TARGET="$MACOSX_DEPLOYMENT_TARGET" \
    SOURCE_DATE_EPOCH="$SOURCE_DATE_EPOCH" RUSTFLAGS="$mac_rustflags" \
    "$rustup" run "$MAC_TOOLCHAIN" cargo build --release --locked -p bammetrics-cli --features pdf --target "$target"
}

mac_universal() {
  local dest="$repo/$work/$(name_of universal-apple-darwin)"
  [ "$dry" = 1 ] || mkdir -p "$dest"
  run lipo -create -output "$dest/bammetrics" \
    "$repo/target/release-mac/aarch64-apple-darwin/release/bammetrics" \
    "$repo/target/release-mac/x86_64-apple-darwin/release/bammetrics"
}

mac_check() {
  local bin="$repo/$work/$(name_of universal-apple-darwin)/bammetrics" got archs
  if [ "$dry" = 1 ]; then
    echo "+ self-check universal-apple-darwin: lipo -archs, --version --json native and arch -x86_64"
    return
  fi
  archs="$(lipo -archs "$bin")"
  [ "$archs" = "x86_64 arm64" ] || fail "universal-apple-darwin has archs '$archs', expected 'x86_64 arm64'"
  got="$("$bin" --version --json)" || fail "universal-apple-darwin arm64 slice does not run"
  version_check "universal-apple-darwin (arm64)" "$got"
  echo "self-check universal-apple-darwin arm64: $got"
  got="$(arch -x86_64 "$bin" --version --json)" || fail "universal-apple-darwin x86_64 slice does not run"
  version_check "universal-apple-darwin (x86_64)" "$got"
  echo "self-check universal-apple-darwin x86_64: $got"
}

for t in x86_64-unknown-linux-musl aarch64-unknown-linux-musl; do
  if has "$t"; then
    start=$SECONDS
    linux_build "$t"
    linux_check "$t"
    echo "release.sh: $t built in $((SECONDS - start)) s" >&2
  fi
done
if has universal-apple-darwin; then
  start=$SECONDS
  mac_build aarch64-apple-darwin
  mac_build x86_64-apple-darwin
  mac_universal
  mac_check
  echo "release.sh: universal-apple-darwin built in $((SECONDS - start)) s" >&2
fi

# --- tarballs, SBOM, sums, formula, build info ------------------------------
tarball() {
  local target="$1" name tar_script
  name="$(name_of "$target")"
  if [ "$dry" = 0 ]; then
    cp README.md CHANGELOG.md LICENSE-MIT LICENSE-APACHE "$work/$name/"
    chmod 755 "$work/$name/bammetrics"
    chmod 644 "$work/$name/README.md" "$work/$name/CHANGELOG.md" "$work/$name/LICENSE-MIT" "$work/$name/LICENSE-APACHE"
  fi
  tar_script="apk add --no-cache tar >/dev/null; cd /out"
  tar_script+="; tar --sort=name --mtime=@$SOURCE_DATE_EPOCH --owner=0 --group=0 --numeric-owner -cf - -C .work $name | gzip -n > $name.tar.gz"
  tar_script+="; chown $(id -u):$(id -g) $name.tar.gz"
  run docker run --rm -v "$repo/$out:/out" "$RUNNER_IMAGE" sh -euc "$tar_script"
}
for t in "${targets[@]}"; do tarball "$t"; done

sbom="bammetrics-$version.cdx.json"
sbom_build() {
  local script
  script="apk add --no-cache musl-dev git >/dev/null; git config --global --add safe.directory '*'"
  script+="; $(pin_tool cargo-cyclonedx "$CARGO_CYCLONEDX_VERSION")"
  script+="; export PATH=/tools/bin:\$PATH"
  script+="; cargo cyclonedx --manifest-path crates/cli/Cargo.toml --format json --spec-version 1.5 --features pdf --override-filename bammetrics-$version"
  script+="; mv crates/cli/bammetrics-$version.json $out/$sbom"
  script+="; rm -f crates/*/bammetrics-$version.json"
  script+="; chown $(id -u):$(id -g) $out/$sbom"
  run docker run --rm \
    -v "$repo:/build" -w /build \
    -v "$git_common:$git_common:ro" \
    -v bammetrics-cargo-registry:/usr/local/cargo/registry \
    -v bammetrics-cargo-tools:/tools \
    -e CARGO_INSTALL_ROOT=/tools -e CARGO_TARGET_DIR=/target \
    "$RUST_ALPINE_IMAGE" sh -euc "$script"
}
sbom_build
# cargo cyclonedx stamps the run's clock and a random serial number; pin both to the commit.
run python3 "$repo/scripts/sbom_pin.py" "$out/$sbom" "$full_commit" "$SOURCE_DATE_EPOCH"

if [ "$dry" = 1 ]; then
  echo "+ write $out/BUILDINFO"
else
  {
    echo "commit $full_commit"
    echo "source_date_epoch $SOURCE_DATE_EPOCH"
    echo "rust_alpine_image $RUST_ALPINE_IMAGE"
    echo "runner_image $RUNNER_IMAGE"
    echo "mac_toolchain $MAC_TOOLCHAIN"
    echo "cargo_auditable $CARGO_AUDITABLE_VERSION"
    echo "cargo_cyclonedx $CARGO_CYCLONEDX_VERSION"
    [ -f "$work/musl-dev.txt" ] && echo "apk $(cat "$work/musl-dev.txt")"
    echo "macosx_deployment_target $MACOSX_DEPLOYMENT_TARGET"
    for t in x86_64-unknown-linux-musl aarch64-unknown-linux-musl; do
      has "$t" && echo "rustc $t $(cat "$work/rustc-$t.txt")"
    done
    has universal-apple-darwin \
      && echo "rustc universal-apple-darwin $("$rustup" run "$MAC_TOOLCHAIN" rustc --version)"
  } > "$out/BUILDINFO"
fi

# No artifact may name the build host: the repository path or $HOME.
leak_check() {
  local scan="$work/scan" needle hit f
  if [ "$dry" = 1 ]; then
    echo "+ leak check: $* for the repository path and \$HOME"
    return
  fi
  rm -rf "$scan"; mkdir -p "$scan"
  for f in "$@"; do
    case "$f" in
      *.tar.gz) mkdir -p "$scan/$(basename "$f")"; tar -xzf "$f" -C "$scan/$(basename "$f")" ;;
      *) cp "$f" "$scan/" ;;
    esac
  done
  for needle in "$repo" "${HOME:-}"; do
    case "$needle" in ""|/) continue ;; esac
    hit="$(grep -rlaF -- "$needle" "$scan" | head -n 1 || true)"
    [ -z "$hit" ] || fail "$needle found in ${hit#"$scan"/}"
  done
  rm -rf "$scan"
}

artifacts=()
for t in "${targets[@]}"; do artifacts+=("$out/$(name_of "$t").tar.gz"); done
artifacts+=("$out/$sbom" "$out/BUILDINFO")
leak_check ${artifacts[@]+"${artifacts[@]}"}

# SHA256SUMS lists what this run made: the tarballs of the unblocked targets and the SBOM.
listed=()
for t in "${targets[@]}"; do listed+=("$(name_of "$t").tar.gz"); done
listed+=("$sbom")
if [ "$dry" = 1 ]; then
  echo "+ (cd $out && shasum -a 256 ${listed[*]} > SHA256SUMS)"
else
  (
    cd "$out"
    # shellcheck disable=SC2046
    shasum -a 256 $(printf '%s\n' "${listed[@]}" | LC_ALL=C sort) > SHA256SUMS
  )
  (cd "$out" && shasum -a 256 -c SHA256SUMS >/dev/null) || fail "SHA256SUMS does not verify"
fi

skip_args=()
for t in x86_64-unknown-linux-musl aarch64-unknown-linux-musl universal-apple-darwin; do
  has "$t" || skip_args+=(--skip "$t")
done
if [ "${#targets[@]}" -gt 0 ]; then
  run python3 "$repo/scripts/brew_formula.py" "$version" "$out/SHA256SUMS" "$out/bammetrics.rb" ${skip_args[@]+"${skip_args[@]}"}
fi

if [ -f "$out/bammetrics.rb" ]; then leak_check "$out/bammetrics.rb"; fi

# --- publish ----------------------------------------------------------------
if [ "$publish" = 1 ]; then
  notes="$work/notes.md"
  if [ "$dry" = 0 ]; then printf '%s\n' "$notes_text" > "$notes"; fi
  files=()
  for t in "${targets[@]}"; do files+=("$out/$(name_of "$t").tar.gz"); done
  files+=("$out/$sbom")
  # A draft first: a failed upload leaves no public release with missing files.
  run gh release create "v$version" --draft --verify-tag --title "BAMmetrics $version" \
    --notes-file "$notes" "${files[@]}" "$out/SHA256SUMS" "$out/BUILDINFO"
  run gh release edit "v$version" --draft=false
fi

[ "$dry" = 1 ] || rm -rf "$work"
echo "release.sh: done, $out" >&2

#!/usr/bin/env python3
"""Write the Homebrew formula bammetrics.rb from a release's SHA256SUMS.

    brew_formula.py VERSION SHA256SUMS OUT [--skip TARGET ...]

Exit 2 when a target that is not skipped has no tarball line in SHA256SUMS.
The formula has no license line while the repository has no licence.
"""
import re
import sys
from pathlib import Path

REPO = "https://github.com/nuin/bammetrics"
MAC = "universal-apple-darwin"
INTEL = "x86_64-unknown-linux-musl"
ARM = "aarch64-unknown-linux-musl"


def die(message):
    print(f"brew_formula.py: {message}", file=sys.stderr)
    sys.exit(2)


def read_sums(path, version):
    sums = {}
    pattern = re.compile(rf"^([0-9a-f]{{64}})  bammetrics-{re.escape(version)}-(.+)\.tar\.gz$")
    for line in Path(path).read_text().splitlines():
        match = pattern.match(line)
        if match:
            sums[match.group(2)] = match.group(1)
    return sums


def block(target, sums, version, indent):
    pad = " " * indent
    url = f"{REPO}/releases/download/v{version}/bammetrics-{version}-{target}.tar.gz"
    return f'{pad}url "{url}"\n{pad}sha256 "{sums[target]}"\n'


def formula(version, sums, skip):
    out = [
        "class Bammetrics < Formula",
        '  desc "Picard BAM QC metrics, reimplemented in Rust"',
        f'  homepage "{REPO}"',
        f'  version "{version}"',
        "",
    ]
    if MAC not in skip:
        # brew style allows url only inside on_arm and on_intel, so both get the universal file.
        out.append("  on_macos do")
        for name in ("on_arm", "on_intel"):
            out += [f"    {name} do", block(MAC, sums, version, 6).rstrip("\n"), "    end"]
        out += ["  end", ""]
    linux = [(n, t) for n, t in (("on_intel", INTEL), ("on_arm", ARM)) if t not in skip]
    if linux:
        out.append("  on_linux do")
        for name, target in linux:
            out += [f"    {name} do", block(target, sums, version, 6).rstrip("\n"), "    end"]
        out += ["  end", ""]
    out += [
        "  def install",
        '    bin.install "bammetrics"',
        "  end",
        "",
        "  test do",
        f'    assert_match "{version}", shell_output("#{{bin}}/bammetrics --version")',
        "  end",
        "end",
        "",
    ]
    return "\n".join(out)


def main(argv):
    skip = set()
    args = []
    it = iter(argv)
    for a in it:
        if a == "--skip":
            skip.add(next(it, ""))
        else:
            args.append(a)
    if len(args) != 3:
        die("usage: brew_formula.py VERSION SHA256SUMS OUT [--skip TARGET ...]")
    version, sums_path, out = args
    sums = read_sums(sums_path, version)
    for target in (MAC, INTEL, ARM):
        if target not in skip and target not in sums:
            die(f"{sums_path} has no tarball line for {target}")
    Path(out).write_text(formula(version, sums, skip))


if __name__ == "__main__":
    main(sys.argv[1:])

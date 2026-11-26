//! Sets `BAMMETRICS_COMMIT` to the short commit the build was made from,
//! for the tool listing. Never fails: without git (a source archive), or
//! when the git checkout around the sources is not this workspace's own (a
//! crate unpacked inside another repository), the variable is left unset
//! and the listing says `null`.

use std::path::Path;
use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let text = String::from_utf8(out.stdout).ok()?;
    let text = text.trim();
    (!text.is_empty()).then(|| text.to_owned())
}

/// True when the repository `git` finds is the workspace's own: its top
/// level is `CARGO_MANIFEST_DIR/../..`.
fn own_checkout() -> bool {
    let Some(top) = git(&["rev-parse", "--show-toplevel"]) else {
        return false;
    };
    let Ok(manifest) = std::env::var("CARGO_MANIFEST_DIR") else {
        return false;
    };
    let root = Path::new(&manifest).join("../..");
    match (Path::new(&top).canonicalize(), root.canonicalize()) {
        (Ok(top), Ok(root)) => top == root,
        _ => false,
    }
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    if !own_checkout() {
        return;
    }
    if let Some(dir) = git(&["rev-parse", "--git-dir"]) {
        let dir = Path::new(&dir);
        println!("cargo:rerun-if-changed={}", dir.join("HEAD").display());
        if let Some(head) = git(&["symbolic-ref", "-q", "HEAD"]) {
            let common = git(&["rev-parse", "--git-common-dir"]);
            let base = common.as_deref().map_or(dir, Path::new);
            println!("cargo:rerun-if-changed={}", base.join(head).display());
        }
    }
    if let Some(commit) = git(&["rev-parse", "--short", "HEAD"]) {
        println!("cargo:rustc-env=BAMMETRICS_COMMIT={commit}");
    }
}

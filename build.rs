/// build.rs
use std::process::Command;

fn main() {
    let output = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    let git_hash = String::from_utf8(output.stdout).unwrap();
    let version = env!("CARGO_PKG_VERSION");
    println!(
        "cargo:rustc-env=CARGO_GIT_HASH={}
        cargo:rustc-env=CARGO_LONG_VERSION={} commit:{}",
        git_hash, version, git_hash
    );
}

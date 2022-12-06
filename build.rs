/// build.rs
use std::process::Command;
fn main() {
    // note: add error checking yourself.
    let output = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    let git_hash = String::from_utf8(output.stdout).unwrap();
    let version = env!("CARGO_PKG_VERSION");
    println!("cargo:rustc-env=GIT_HASH=v{} commit:{}'", version, git_hash);
}

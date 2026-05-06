use std::path::{Path, PathBuf};
use std::process::Command;

pub fn ft() -> PathBuf {
    env!("CARGO_BIN_EXE_ft").into()
}

pub fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// Run ft with the given args; return stdout as a String. Panics on non-zero exit.
pub fn run(args: &[&str]) -> String {
    let out = Command::new(ft())
        .args(args)
        .output()
        .expect("failed to spawn ft");
    assert!(
        out.status.success(),
        "ft exited {}\nstderr: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("non-UTF8 stdout")
}

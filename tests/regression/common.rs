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

/// Select named columns from TSV output by header name, in the order given.
/// Snapshots only the specified columns so that adding new columns to a
/// command's output doesn't count as a regression.
pub fn select_tsv_cols(tsv: &str, cols: &[&str]) -> String {
    let mut lines = tsv.lines();
    let headers: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let indices: Vec<usize> = cols
        .iter()
        .map(|c| {
            headers
                .iter()
                .position(|h| h == c)
                .unwrap_or_else(|| panic!("column {c:?} not found; headers: {headers:?}"))
        })
        .collect();
    let mut out = cols.join("\t") + "\n";
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        out += &indices
            .iter()
            .map(|&i| fields.get(i).copied().unwrap_or(""))
            .collect::<Vec<_>>()
            .join("\t");
        out += "\n";
    }
    out
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

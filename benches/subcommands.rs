//! End-to-end criterion benches for the read-heavy fibertools subcommands
//! (`pileup`, `fire`, `track-decorators`, `extract`). Each emits per-record
//! annotations (m6a / cpg / msp / nuc / fire), so wall-clock is dominated by
//! annotation traversal cost.
//!
//! Each bench shells out to the freshly-built `ft` binary so we measure what
//! the user actually runs. Inputs come from the public FIRE test-data CRAM;
//! fetch via `bash benches/fetch-data.sh` before running.

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use tempfile::TempDir;

const FT_BIN: &str = env!("CARGO_BIN_EXE_ft");

fn data_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("benches/data/fire-test-data")
}

fn cram() -> PathBuf {
    let p = data_root().join("test.cram");
    assert!(
        p.exists(),
        "bench data missing at {} \u{2014} run `bash benches/fetch-data.sh` first",
        p.display()
    );
    p
}

fn run_ft(args: &[&str], cwd: &Path) {
    let status = Command::new(FT_BIN)
        .env("REF_PATH", data_root().join("test.fa.gz"))
        .args(args)
        .current_dir(cwd)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("failed to spawn ft");
    assert!(
        status.success(),
        "ft exited non-zero ({:?}) with args: {:?}",
        status.code(),
        args
    );
}

/// Run `ft <args>` once per iteration in a fresh tempdir (so output files don't
/// accumulate). Output paths in `args` are relative — they land in the tempdir.
fn bench_cmd(c: &mut Criterion, name: &str, args: &[&str]) {
    let mut g = c.benchmark_group(name);
    g.sample_size(10);
    g.bench_function("default", |b| {
        b.iter_batched(
            || TempDir::new().expect("tempdir creation failed"),
            |tmp| run_ft(args, tmp.path()),
            BatchSize::PerIteration,
        );
    });
    g.finish();
}

fn bench_pileup(c: &mut Criterion) {
    let cram = cram();
    let cram = cram.to_str().expect("cram path is not utf-8");
    bench_cmd(
        c,
        "pileup",
        &["pileup", "--m6a", "--cpg", "-o", "pileup.bed.gz", cram],
    );
}

fn bench_fire(c: &mut Criterion) {
    let cram = cram();
    let cram = cram.to_str().expect("cram path is not utf-8");
    bench_cmd(c, "fire", &["fire", cram, "fire.bam"]);
}

fn bench_decorator(c: &mut Criterion) {
    let cram = cram();
    let cram = cram.to_str().expect("cram path is not utf-8");
    bench_cmd(
        c,
        "decorator",
        &[
            "track-decorators",
            "--bed12",
            "decorator.bed12",
            "--decorator",
            "decorator.dec.bed",
            cram,
        ],
    );
}

fn bench_extract(c: &mut Criterion) {
    let cram = cram();
    let cram = cram.to_str().expect("cram path is not utf-8");
    bench_cmd(c, "extract", &["extract", "--all", "extract.tsv.gz", cram]);
}

criterion_group!(
    benches,
    bench_pileup,
    bench_fire,
    bench_decorator,
    bench_extract
);
criterion_main!(benches);

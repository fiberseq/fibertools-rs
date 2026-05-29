//! End-to-end criterion benches for the read-heavy fibertools subcommands
//! (`pileup`, `fire`, `track-decorators`, `extract`). Each one iterates or
//! emits per-record annotations (m6a / cpg / msp / nuc / fire) on the input
//! BAM/CRAM, so wall-clock here is dominated by annotation traversal cost.
//!
//! Each bench shells out to the freshly-built `ft` binary
//! (`env!("CARGO_BIN_EXE_ft")`) so we measure what the user actually runs.
//! Inputs come from the public FIRE test-data CRAM; fetch via
//! `bash benches/fetch-data.sh` before running.

use criterion::measurement::WallTime;
use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkGroup, Criterion,
};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Duration;
use tempfile::TempDir;

const FT_BIN: &str = env!("CARGO_BIN_EXE_ft");

fn data_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("benches/data/fire-test-data")
}

fn cram() -> PathBuf {
    data_root().join("test.cram")
}

fn ref_fa() -> PathBuf {
    data_root().join("test.fa.gz")
}

fn assert_data_present() {
    let p = cram();
    assert!(
        p.exists(),
        "bench data missing at {} \u{2014} run `bash benches/fetch-data.sh` first",
        p.display()
    );
}

fn configure<'a>(c: &'a mut Criterion, name: &str) -> BenchmarkGroup<'a, WallTime> {
    let mut g = c.benchmark_group(name);
    g.sample_size(10);
    g.measurement_time(Duration::from_secs(60));
    g
}

fn run_ft(args: &[&str], cwd: &Path) {
    let status = Command::new(FT_BIN)
        .env("REF_PATH", ref_fa())
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

fn bench_pileup(c: &mut Criterion) {
    assert_data_present();
    let mut g = configure(c, "pileup");
    // pileup over the whole FIRE test-data CRAM exceeds the default
    // measurement budget; give it more headroom rather than scoping with --rgn
    // so we still measure real-world wall-clock.
    g.measurement_time(Duration::from_secs(180));
    let cram_path = cram();
    let cram_str = cram_path.to_str().expect("cram path is not utf-8");

    g.bench_function("default", |b| {
        b.iter_batched(
            || TempDir::new().expect("tempdir creation failed"),
            |tmp| {
                let out = tmp.path().join("pileup.bed.gz");
                let out_str = out.to_str().expect("out path is not utf-8");
                run_ft(
                    &["pileup", "--m6a", "--cpg", "-o", out_str, cram_str],
                    tmp.path(),
                );
            },
            BatchSize::PerIteration,
        );
    });

    g.finish();
}

fn bench_fire(c: &mut Criterion) {
    assert_data_present();
    let mut g = configure(c, "fire");
    // fire over the whole FIRE test-data CRAM exceeds the default
    // measurement budget; raise it so we still measure real-world wall-clock.
    g.measurement_time(Duration::from_secs(120));
    let cram_path = cram();
    let cram_str = cram_path.to_str().expect("cram path is not utf-8");

    g.bench_function("default", |b| {
        b.iter_batched(
            || TempDir::new().expect("tempdir creation failed"),
            |tmp| {
                let out = tmp.path().join("fire.bam");
                let out_str = out.to_str().expect("out path is not utf-8");
                run_ft(&["fire", cram_str, out_str], tmp.path());
            },
            BatchSize::PerIteration,
        );
    });

    g.finish();
}

fn bench_decorator(c: &mut Criterion) {
    assert_data_present();
    let mut g = configure(c, "decorator");
    let cram_path = cram();
    let cram_str = cram_path.to_str().expect("cram path is not utf-8");

    g.bench_function("default", |b| {
        b.iter_batched(
            || TempDir::new().expect("tempdir creation failed"),
            |tmp| {
                let bed12 = tmp.path().join("decorator.bed12");
                let dec = tmp.path().join("decorator.dec.bed");
                let bed12_str = bed12.to_str().expect("bed12 path is not utf-8");
                let dec_str = dec.to_str().expect("dec path is not utf-8");
                run_ft(
                    &[
                        "track-decorators",
                        "--bed12",
                        bed12_str,
                        "--decorator",
                        dec_str,
                        cram_str,
                    ],
                    tmp.path(),
                );
            },
            BatchSize::PerIteration,
        );
    });

    g.finish();
}

fn bench_extract(c: &mut Criterion) {
    assert_data_present();
    let mut g = configure(c, "extract");
    let cram_path = cram();
    let cram_str = cram_path.to_str().expect("cram path is not utf-8");

    g.bench_function("all", |b| {
        b.iter_batched(
            || TempDir::new().expect("tempdir creation failed"),
            |tmp| {
                let out = tmp.path().join("extract.tsv.gz");
                let out_str = out.to_str().expect("out path is not utf-8");
                run_ft(&["extract", "--all", out_str, cram_str], tmp.path());
            },
            BatchSize::PerIteration,
        );
    });

    g.finish();
}

criterion_group!(
    benches,
    bench_pileup,
    bench_fire,
    bench_decorator,
    bench_extract
);
criterion_main!(benches);

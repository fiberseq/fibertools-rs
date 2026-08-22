use super::common::{fixture, run};

/// Pin ft call-peaks output so refactors of the shared peak caller
/// (call_peaks_for_chrom is also union-peaks' engine) cannot drift it.
#[test]
fn call_peaks_ctcf_snapshot() {
    let out = run(&[
        "call-peaks",
        fixture("ctcf.bam").to_str().unwrap(),
        "--min-fire-frac",
        "0.5",
    ]);
    insta::assert_snapshot!(out);
}

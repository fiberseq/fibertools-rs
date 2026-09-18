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

/// `--haps` must fill the H1/H2 columns from HP tags (issue #140: the flag was
/// parsed but never reached the pileup, so every haplotype column was zero).
#[test]
fn call_peaks_haps_fills_haplotype_columns() {
    let out = run(&[
        "call-peaks",
        fixture("NAPA.bam").to_str().unwrap(),
        "--haps",
        "--min-fire-frac",
        "0.5",
    ]);
    let peak = out.lines().find(|l| !l.starts_with('#')).expect("one peak");
    let cols: Vec<&str> = peak.split('\t').collect();
    // coverage, coverage_H1, coverage_H2
    assert_eq!(cols[5], "95");
    assert_eq!(cols[10], "45");
    assert_eq!(cols[15], "14");
}

use super::common::{fixture, run};

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

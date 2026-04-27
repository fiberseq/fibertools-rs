//! Unit tests for the `call_peaks` module's pure helpers and FDR machinery.
//!
//! These tests exercise public surface only — `lookup_fdr`,
//! `reciprocal_overlap_raw`, `IncrementalFdrBuilder`, and the FDR-table
//! read/write helpers. The geometry-only `reciprocal_overlap_raw` and the
//! score-only `lookup_fdr` were extracted from `Peak` methods so they can
//! be tested without constructing a `FiberseqPileup`.

use fibertools_rs::subcommands::call_peaks::{
    lookup_fdr, read_fdr_table, reciprocal_overlap_raw, write_fdr_table, FdrEntry,
    IncrementalFdrBuilder, PileupRecord,
};
use std::io::Write;
use tempfile::NamedTempFile;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

fn fdr_entry(threshold: f64, fdr: f64, shuffled_bp: f64, real_bp: f64) -> FdrEntry {
    FdrEntry {
        threshold,
        fdr,
        shuffled_bp,
        real_bp,
    }
}

fn pileup_record(start: u64, end: u64, score: f64) -> PileupRecord {
    PileupRecord {
        start,
        end,
        coverage: 1,
        fire_coverage: 1,
        score,
    }
}

// ---------------------------------------------------------------------------
// lookup_fdr
// ---------------------------------------------------------------------------

#[test]
fn lookup_fdr_empty_table_returns_one() {
    let table: Vec<FdrEntry> = vec![];
    assert_eq!(lookup_fdr(0.0, &table), 1.0);
    assert_eq!(lookup_fdr(99.0, &table), 1.0);
    assert_eq!(lookup_fdr(-5.0, &table), 1.0);
}

#[test]
fn lookup_fdr_score_below_all_thresholds_returns_first_entry_fdr() {
    let table = vec![
        fdr_entry(1.0, 0.5, 0.0, 0.0),
        fdr_entry(5.0, 0.1, 0.0, 0.0),
        fdr_entry(10.0, 0.01, 0.0, 0.0),
    ];
    // Score below the smallest threshold returns the first entry's FDR.
    assert!(approx_eq(lookup_fdr(0.5, &table), 0.5));
    assert!(approx_eq(lookup_fdr(-100.0, &table), 0.5));
}

#[test]
fn lookup_fdr_score_above_all_thresholds_returns_last_entry_fdr() {
    let table = vec![
        fdr_entry(1.0, 0.5, 0.0, 0.0),
        fdr_entry(5.0, 0.1, 0.0, 0.0),
        fdr_entry(10.0, 0.01, 0.0, 0.0),
    ];
    assert!(approx_eq(lookup_fdr(100.0, &table), 0.01));
}

#[test]
fn lookup_fdr_exact_threshold_match_returns_that_fdr() {
    let table = vec![
        fdr_entry(1.0, 0.5, 0.0, 0.0),
        fdr_entry(5.0, 0.1, 0.0, 0.0),
        fdr_entry(10.0, 0.01, 0.0, 0.0),
    ];
    assert!(approx_eq(lookup_fdr(1.0, &table), 0.5));
    assert!(approx_eq(lookup_fdr(5.0, &table), 0.1));
    assert!(approx_eq(lookup_fdr(10.0, &table), 0.01));
}

#[test]
fn lookup_fdr_score_between_thresholds_returns_largest_le_score() {
    let table = vec![
        fdr_entry(1.0, 0.5, 0.0, 0.0),
        fdr_entry(5.0, 0.1, 0.0, 0.0),
        fdr_entry(10.0, 0.01, 0.0, 0.0),
    ];
    // Between 1.0 and 5.0 → uses the 1.0 entry.
    assert!(approx_eq(lookup_fdr(3.0, &table), 0.5));
    assert!(approx_eq(lookup_fdr(4.99, &table), 0.5));
    // Between 5.0 and 10.0 → uses the 5.0 entry.
    assert!(approx_eq(lookup_fdr(7.0, &table), 0.1));
    assert!(approx_eq(lookup_fdr(9.99, &table), 0.1));
}

#[test]
fn lookup_fdr_is_monotone_nonincreasing_for_decreasing_fdr_table() {
    // Real FDR tables are sorted ascending by threshold and FDR is generally
    // non-increasing as threshold rises. lookup_fdr should preserve this:
    // higher score → equal-or-lower FDR.
    let table = vec![
        fdr_entry(0.0, 0.9, 0.0, 0.0),
        fdr_entry(1.0, 0.5, 0.0, 0.0),
        fdr_entry(2.5, 0.2, 0.0, 0.0),
        fdr_entry(5.0, 0.05, 0.0, 0.0),
        fdr_entry(10.0, 0.001, 0.0, 0.0),
    ];
    let mut prev = f64::INFINITY;
    for s in [-1.0, 0.0, 0.5, 1.0, 2.0, 2.5, 3.7, 5.0, 8.0, 10.0, 50.0] {
        let f = lookup_fdr(s as f32, &table);
        assert!(
            f <= prev + 1e-12,
            "non-monotonic at score {s}: {f} > {prev}"
        );
        prev = f;
    }
}

// ---------------------------------------------------------------------------
// reciprocal_overlap_raw
// ---------------------------------------------------------------------------

#[test]
fn reciprocal_overlap_different_chromosomes_is_zero() {
    assert_eq!(reciprocal_overlap_raw("chr1", 0, 100, "chr2", 0, 100), 0.0);
}

#[test]
fn reciprocal_overlap_disjoint_intervals_is_zero() {
    assert_eq!(
        reciprocal_overlap_raw("chr1", 0, 100, "chr1", 200, 300),
        0.0
    );
}

#[test]
fn reciprocal_overlap_touching_boundaries_is_zero() {
    // BED-style half-open: end == other.start means no overlapping bases.
    assert_eq!(
        reciprocal_overlap_raw("chr1", 0, 100, "chr1", 100, 200),
        0.0
    );
}

#[test]
fn reciprocal_overlap_identical_intervals_is_one() {
    assert!(approx_eq(
        reciprocal_overlap_raw("chr1", 100, 200, "chr1", 100, 200),
        1.0
    ));
}

#[test]
fn reciprocal_overlap_full_containment_returns_smaller_over_larger() {
    // [100, 200] contains [120, 180]: overlap = 60, smaller_len = 60, larger_len = 100.
    // Reciprocal overlap = min(60/100, 60/60) = 0.6.
    assert!(approx_eq(
        reciprocal_overlap_raw("chr1", 100, 200, "chr1", 120, 180),
        0.6
    ));
}

#[test]
fn reciprocal_overlap_partial_overlap_returns_min_fraction() {
    // [0, 100] and [50, 200]: overlap = 50.
    // a_frac = 50/100 = 0.5, b_frac = 50/150 ≈ 0.333; min ≈ 0.333.
    let v = reciprocal_overlap_raw("chr1", 0, 100, "chr1", 50, 200);
    assert!(approx_eq(v, 50.0 / 150.0));
}

#[test]
fn reciprocal_overlap_is_symmetric() {
    let a = reciprocal_overlap_raw("chr1", 0, 100, "chr1", 50, 200);
    let b = reciprocal_overlap_raw("chr1", 50, 200, "chr1", 0, 100);
    assert!(approx_eq(a, b));
}

#[test]
fn reciprocal_overlap_fifty_percent_each_side() {
    // [0, 100] and [50, 150]: overlap = 50, both lens = 100.
    assert!(approx_eq(
        reciprocal_overlap_raw("chr1", 0, 100, "chr1", 50, 150),
        0.5
    ));
}

// ---------------------------------------------------------------------------
// IncrementalFdrBuilder
// ---------------------------------------------------------------------------

#[test]
fn fdr_builder_single_chromosome_yields_expected_table() {
    // One real record (score=1.0, bp=200) and one shuffled record (score=0.5, bp=100).
    // After fdr_from_fire_scores in descending-by-score order:
    //   First we push when score changes from 1.0 → 0.5 with cur_r=200, cur_v=0.
    //   Then we accumulate cur_v=100, then push the trailing sentinel
    //   (1.0 shuffled, 1.0 real, threshold=-1.0).
    //   FDRs: 0/200 = 0.0 and 1/1 = 1.0.
    let real = vec![pileup_record(0, 200, 1.0)];
    let shuffled = vec![pileup_record(0, 100, 0.5)];

    let mut builder = IncrementalFdrBuilder::new();
    builder.add_chromosome_data(&real, &shuffled);
    let table = builder.build(0.5).expect("FDR build should succeed");

    // Expect two entries sorted ascending by threshold: -1.00 and 1.00.
    assert_eq!(table.len(), 2, "table = {:?}", table);
    assert!(approx_eq(table[0].threshold, -1.0));
    assert!(approx_eq(table[0].fdr, 1.0));
    assert!(approx_eq(table[1].threshold, 1.0));
    assert!(approx_eq(table[1].fdr, 0.0));
}

#[test]
fn fdr_builder_bails_when_no_threshold_below_max_fdr() {
    // Only shuffled data → every FDR is 1.0, so any max_fdr < 1.0 must error.
    let real: Vec<PileupRecord> = vec![];
    let shuffled = vec![pileup_record(0, 100, 1.0)];

    let mut builder = IncrementalFdrBuilder::new();
    builder.add_chromosome_data(&real, &shuffled);
    let err = builder
        .build(0.5)
        .expect_err("should bail with no real data");
    let msg = format!("{err}");
    assert!(
        msg.contains("FDR"),
        "expected error message to mention FDR, got: {msg}"
    );
}

#[test]
fn fdr_builder_multi_chromosome_additivity_matches_single_call() {
    // Splitting the same input across two add_chromosome_data calls must
    // produce the same final table as feeding it all at once. The builder
    // aggregates by score key so per-chromosome boundaries should not matter.
    let real_a = vec![pileup_record(0, 200, 1.0), pileup_record(0, 50, 2.0)];
    let shuffled_a = vec![pileup_record(0, 100, 0.5)];
    let real_b = vec![pileup_record(0, 75, 1.0)];
    let shuffled_b = vec![pileup_record(0, 25, 0.5), pileup_record(0, 60, 1.0)];

    let mut combined = IncrementalFdrBuilder::new();
    let real_all: Vec<_> = real_a.iter().chain(real_b.iter()).cloned().collect();
    let shuffled_all: Vec<_> = shuffled_a
        .iter()
        .chain(shuffled_b.iter())
        .cloned()
        .collect();
    combined.add_chromosome_data(&real_all, &shuffled_all);
    let table_combined = combined.build(0.99).expect("build should succeed");

    let mut split = IncrementalFdrBuilder::new();
    split.add_chromosome_data(&real_a, &shuffled_a);
    split.add_chromosome_data(&real_b, &shuffled_b);
    let table_split = split.build(0.99).expect("build should succeed");

    assert_eq!(
        table_combined.len(),
        table_split.len(),
        "combined={:?}\nsplit={:?}",
        table_combined,
        table_split
    );
    for (c, s) in table_combined.iter().zip(table_split.iter()) {
        assert!(approx_eq(c.threshold, s.threshold));
        assert!(approx_eq(c.fdr, s.fdr));
        assert!(approx_eq(c.shuffled_bp, s.shuffled_bp));
        assert!(approx_eq(c.real_bp, s.real_bp));
    }
}

#[test]
fn fdr_builder_output_is_sorted_ascending_by_threshold() {
    let real = vec![
        pileup_record(0, 200, 5.0),
        pileup_record(0, 100, 2.0),
        pileup_record(0, 50, 1.0),
    ];
    let shuffled = vec![pileup_record(0, 100, 4.0), pileup_record(0, 100, 0.5)];

    let mut builder = IncrementalFdrBuilder::new();
    builder.add_chromosome_data(&real, &shuffled);
    let table = builder.build(0.99).expect("build should succeed");

    let mut prev = f64::NEG_INFINITY;
    for entry in &table {
        assert!(
            entry.threshold > prev || approx_eq(entry.threshold, prev),
            "thresholds not ascending: {} after {}",
            entry.threshold,
            prev
        );
        prev = entry.threshold;
    }
}

// ---------------------------------------------------------------------------
// read_fdr_table / write_fdr_table round trip
// ---------------------------------------------------------------------------

#[test]
fn fdr_table_round_trip_via_tempfile() {
    let original = vec![
        fdr_entry(-1.0, 1.0, 1.0, 1.0),
        fdr_entry(0.5, 0.25, 100.0, 400.0),
        fdr_entry(1.0, 0.05, 50.0, 1000.0),
    ];

    let tmp = NamedTempFile::new().expect("tempfile");
    let path = tmp.path().to_str().unwrap().to_string();
    write_fdr_table(&original, &path).expect("write");
    let round_tripped = read_fdr_table(&path).expect("read");

    assert_eq!(round_tripped.len(), original.len());
    for (a, b) in round_tripped.iter().zip(original.iter()) {
        // write_fdr_table formats threshold/fdr/bp at 2/6/0 decimal places.
        assert!(approx_eq(a.threshold, b.threshold));
        assert!((a.fdr - b.fdr).abs() < 1e-6);
        assert!(approx_eq(a.shuffled_bp, b.shuffled_bp));
        assert!(approx_eq(a.real_bp, b.real_bp));
    }
}

#[test]
fn read_fdr_table_rejects_wrong_column_count() {
    let mut tmp = NamedTempFile::new().expect("tempfile");
    writeln!(tmp, "threshold\tFDR\tshuffled_bp\treal_bp").unwrap();
    writeln!(tmp, "1.0\t0.5\t10").unwrap(); // only 3 columns
    tmp.flush().unwrap();
    let path = tmp.path().to_str().unwrap();

    let err = read_fdr_table(path).expect_err("should reject malformed line");
    let msg = format!("{err}");
    assert!(
        msg.contains("4 columns") || msg.contains("Invalid FDR table format"),
        "expected column-count error, got: {msg}"
    );
}

#[test]
fn read_fdr_table_returns_entries_sorted_ascending_by_threshold() {
    // write a table with entries in non-ascending order; read should return them sorted.
    let mut tmp = NamedTempFile::new().expect("tempfile");
    writeln!(tmp, "threshold\tFDR\tshuffled_bp\treal_bp").unwrap();
    writeln!(tmp, "5.00\t0.010000\t10\t1000").unwrap();
    writeln!(tmp, "1.00\t0.500000\t100\t200").unwrap();
    writeln!(tmp, "-1.00\t1.000000\t1\t1").unwrap();
    tmp.flush().unwrap();
    let path = tmp.path().to_str().unwrap();

    let table = read_fdr_table(path).expect("read");
    assert_eq!(table.len(), 3);
    assert!(table[0].threshold < table[1].threshold);
    assert!(table[1].threshold < table[2].threshold);
}

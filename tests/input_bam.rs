use fibertools_rs::utils::input_bam::*;

use fibertools_rs::fiber::FiberseqData;
use fibertools_rs::utils::basemods::M6A_TYPE;
use molecular_annotation::{Encoding, MolecularAnnotations, QualitySpec, Strand};

/// Build a minimal `FiberseqData` with the given msp lengths and m6A
/// count, then derive its fiberseq_callable state at the given minimums,
/// the same path the stream drop consults through `is_callable`.
fn make_fsd(msp_lengths: &[i64], m6a_count: usize, filters: &FiberFilters) -> FiberseqData {
    let mut annotations = MolecularAnnotations::new(1000);
    if !msp_lengths.is_empty() {
        let t = annotations.add_annotation_type("msp", QualitySpec::none(), Encoding::Ma);
        let mut start = 0u32;
        for &len in msp_lengths {
            t.add(start, len as u32, Strand::Forward, vec![], None);
            start += len as u32 + 1;
        }
    }
    if m6a_count > 0 {
        let t = annotations.add_annotation_type(
            M6A_TYPE,
            "Q".parse().expect("Q parses"),
            Encoding::mm_ml(),
        );
        for i in 0..m6a_count {
            t.add(i as u32, 1, Strand::Forward, vec![0], None);
        }
    }
    let (min_msp, min_ave) = filters.callable_minimums();
    fibertools_rs::utils::ma_io::derive_fiberseq_callable(&mut annotations, min_msp, min_ave);
    FiberseqData {
        record: rust_htslib::bam::Record::new(),
        annotations,
        ec: 0.0,
        target_name: ".".to_string(),
        rg: ".".to_string(),
        center_position: None,
    }
}

#[test]
fn minimums_never_drop_without_the_flag() {
    let mut f = FiberFilters::default();
    f.min_msp = Some(100);
    assert!(!f.drop_uncallable_fibers, "the minimums only mark callable");
}

#[test]
fn callable_state_matches_the_resolved_minimums() {
    // Defaults (10 MSPs / mean 10): the classic cases.
    let f = FiberFilters::default();
    let ten_long: Vec<i64> = vec![100; 10];
    assert!(make_fsd(&ten_long, 5, &f).is_callable());
    let nine_long: Vec<i64> = vec![100; 9];
    assert!(!make_fsd(&nine_long, 5, &f).is_callable(), "too few MSPs");
    let ten_short: Vec<i64> = vec![5; 10];
    assert!(!make_fsd(&ten_short, 5, &f).is_callable(), "mean too small");
    // No decodable m6A type: judged from the MSPs alone. The nucleosome
    // caller emits no MSP for a read without m6A, so a surviving MSP is
    // itself proof of m6A at calling time; derivation never inspects the
    // m6A type (which cannot decode on SEQ-less records).
    assert!(make_fsd(&ten_long, 0, &f).is_callable());
}

#[test]
fn msp_presence_is_required_even_at_zero_minimums() {
    // The fixed gate cannot be tuned away: even at 0/0 an MSP-less
    // fiber is uncallable. An m6A-less fiber has no MSPs by
    // construction (the caller emits none without m6A), which is why
    // --skip-no-m6a is gone.
    let mut f = FiberFilters::default();
    f.min_msp = Some(0);
    f.min_ave_msp_size = Some(0);
    assert!(!make_fsd(&[], 5, &f).is_callable(), "no MSPs");
    assert!(
        make_fsd(&[100], 1, &f).is_callable(),
        "one MSP passes at 0/0"
    );
}

#[test]
fn custom_minimums_change_the_callable_state() {
    let mut f = FiberFilters::default();
    f.min_msp = Some(5);
    assert!(f.callable_minimums_are_custom());
    assert_eq!(f.callable_minimums(), (5, 10));
    let five_long: Vec<i64> = vec![100; 5];
    assert!(make_fsd(&five_long, 5, &f).is_callable());
    assert!(!make_fsd(&five_long, 5, &FiberFilters::default()).is_callable());
}

use bamlift::*;

///
/// TESTS
///
#[test]
fn test_lift_reference_positions() {
    let aligned_block_pairs = vec![
        ([0, 4], [40, 44]),
        ([7, 8], [50, 51]),
        ([9, 20], [59, 70]),
        ([20, 21], [70, 71]),
    ];
    let positions = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
    let expected_result = vec![41, 42, 43, 44, 44, 50, 50, 51, 59];
    let result = lift_reference_positions(&aligned_block_pairs, &positions);
    assert_eq!(result, expected_result);
}

///```
/// use rust_htslib::{bam, bam::Read};
/// use bamlift::*;
/// use log;
/// use env_logger::{Builder, Target};;
/// Builder::new().target(Target::Stderr).filter(None, log::LevelFilter::Debug).init();
/// let mut bam = bam::Reader::from_path(&"../tests/data/all.bam").unwrap();
/// for record in bam.records() {
///     let record = record.unwrap();
///     let seq_len = i64::try_from(record.seq_len()).unwrap();
///     let positions: Vec<i64> = (0..seq_len).collect();
///     closest_reference_positions(&record, &positions);
/// }
///```
fn _fake() {}

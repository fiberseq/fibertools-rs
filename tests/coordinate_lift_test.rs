use molecular_annotation::liftover::AlignedBlocks;
use rust_htslib::bam;
use rust_htslib::bam::Read;

/// Exercise `AlignedBlocks::lift_to_reference` against a real fiberseq
/// record. Forward query range 3525-3628 (a 103bp nucleosome) should
/// lift to a reference range of the same length when the range sits
/// entirely inside an aligned block.
#[test]
fn test_coordinate_lift() -> anyhow::Result<()> {
    let mut bam = bam::Reader::from_path("tests/data/nuc_example.bam")?;

    let record = bam
        .records()
        .next()
        .expect("nuc_example.bam should have at least one record")?;

    let blocks = AlignedBlocks::from_record(&record);
    let (ref_start, ref_end) = blocks.lift_to_reference(3525, 3628);

    let ref_start = ref_start.expect("3525 should lift");
    let ref_end = ref_end.expect("3628 should lift");
    let ref_length = ref_end - ref_start;

    assert_eq!(
        ref_length, 103,
        "lift should preserve length when the range is inside an aligned block (got {ref_start}-{ref_end})"
    );

    Ok(())
}

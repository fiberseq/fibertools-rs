use molecular_annotation::liftover::AlignedBlocks;

/// Regression for GitHub issue #90.
///
/// Setup:
///              0         1         2
///              01234567890123456789012
/// Read MSP:        | MSP  |
/// Read         ===========-------=====
/// Genome       =======================
/// Genome MSP:      |     MSP     |
///
/// Aligned blocks:
///   Block 1: query [0, 11)  -> ref [0, 11)
///   Block 2: query [11, 16) -> ref [18, 23)  (5bp gap on ref between blocks)
///
/// MSP query range to lift: [4, 11). The end coordinate sits exactly on
/// the boundary between block 1 and block 2.
///
/// The old `lift_reference_positions` (closest-snap) mapped query
/// position 11 to ref 18 (start of block 2), yielding MSP length
/// 18 - 4 = 14 in reference space.
///
/// The spec library's `AlignedBlocks::lift_to_reference` (inward-snap)
/// correctly maps [4, 11) to ref [4, 11) — the requested range is
/// entirely inside block 1, and the half-open end snaps to the last
/// included base + 1 within the block.
#[test]
fn test_alignment_boundary_bug_fixed_by_inward_snap() {
    let blocks = AlignedBlocks::new(
        vec![([0, 11], [0, 11]), ([11, 16], [18, 23])],
        16, // query_len
    );

    let (ref_start, ref_end) = blocks.lift_to_reference(4, 11);
    assert_eq!(
        ref_start,
        Some(4),
        "start of MSP should lift to ref position 4"
    );
    assert_eq!(
        ref_end,
        Some(11),
        "end of MSP at block boundary should lift to ref 11 (end of block 1), \
         not ref 18 (start of block 2)"
    );

    let length = ref_end.unwrap() - ref_start.unwrap();
    assert_eq!(length, 7, "MSP length in reference space");
}

use fibertools_rs::utils::bamlift::lift_reference_positions;

#[test]
fn test_alignment_boundary_bug() {
    // Create alignment blocks that reproduce the issue described in GitHub issue #90
    // Based on the illustration:
    //              0         1         2
    //              01234567890123456789012
    // Read MSP:        | MSP  |
    // Read         ===========-------=====
    // Genome       =======================
    // Genome MSP:      |     MSP     |

    // Alignment blocks:
    // Block 1: query [0, 11) -> reference [0, 11)  (positions 0-10 aligned)
    // Block 2: query [11, 16) -> reference [18, 23) (positions 11-15 -> 18-22, gap 11-17)
    let aligned_block_pairs = vec![
        ([0, 11], [0, 11]),   // First aligned block
        ([11, 16], [18, 23]), // Second aligned block (with gap)
    ];

    // MSP/nucleosome ends exactly at position 11 (end of first alignment block)
    let query_positions = vec![4, 11]; // MSP from position 4 to 11

    println!("Alignment blocks: {:?}", aligned_block_pairs);
    println!("Query positions to lift: {:?}", query_positions);

    // Lift the positions
    let result = lift_reference_positions(&aligned_block_pairs, &query_positions).unwrap();
    println!("Lifted reference positions: {:?}", result);

    // The bug: position 11 should map to reference position 11 (end of first block)
    // but the current code maps it to position 18 (start of second block)
    assert_eq!(result[0], Some(4)); // Start position should be correct

    // This is the problematic case - position 11 at block boundary
    println!("Position 11 maps to reference: {:?}", result[1]);

    // According to the issue, this currently maps to 18 but should map to 11
    // The MSP length becomes 18-4=14 instead of 11-4=7
    if let Some(end_pos) = result[1] {
        let msp_length = end_pos - result[0].unwrap();
        println!("MSP length in reference: {}", msp_length);
        println!("Expected length: 7, actual length: {}", msp_length);

        // The bug manifests as the end position being mapped to the start of the next block
        if end_pos == 18 {
            println!("BUG CONFIRMED: End position {} mapped to start of next block instead of end of current block", end_pos);
        } else if end_pos == 11 {
            println!("FIXED: End position correctly mapped to {}", end_pos);
        }
    }
}

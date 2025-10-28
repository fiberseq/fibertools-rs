use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::Read;

#[test]
fn test_coordinate_lift() -> anyhow::Result<()> {
    let mut bam = bam::Reader::from_path("tests/data/nuc_example.bam")?;

    for result in bam.records() {
        let record = result?;

        // Print basic record info
        println!("Read: {}", String::from_utf8_lossy(record.qname()));
        println!(
            "Reference: {}, Position: {}, End: {}",
            record.tid(),
            record.reference_start(),
            record.reference_end()
        );
        println!("Is reverse: {}", record.is_reverse());
        println!("Sequence length: {}", record.seq_len());

        // Get aligned block pairs for understanding the alignment
        let aligned_blocks: Vec<_> = record.aligned_block_pairs().collect();
        println!("Total aligned blocks: {}", aligned_blocks.len());
        println!(
            "First 5 aligned blocks: {:?}",
            &aligned_blocks[..std::cmp::min(5, aligned_blocks.len())]
        );

        // Look for alignment blocks around position 3525-3628
        println!("Looking for blocks around nucleosome position 3525-3628:");
        for (i, ([q_st, q_en], [r_st, r_en])) in aligned_blocks.iter().enumerate() {
            if (*q_st <= 3628 && *q_en >= 3525) || (*q_st >= 3520 && *q_st <= 3630) {
                println!(
                    "  Block {}: query=[{}, {}) -> ref=[{}, {})",
                    i, q_st, q_en, r_st, r_en
                );
            }
        }

        // Test the problematic nucleosome coordinates
        // Forward position 3525, length 103 -> ends at 3628
        let starts = vec![3525];
        let ends = vec![3628];

        println!("\nTesting nucleosome at forward pos 3525-3628:");

        // Use the lift_query_range function from bamlift
        match fibertools_rs::utils::bamlift::lift_query_range(&record, &starts, &ends) {
            Ok((ref_starts, ref_ends, ref_lengths)) => {
                if let (Some(ref_start), Some(ref_end), Some(ref_length)) =
                    (ref_starts[0], ref_ends[0], ref_lengths[0])
                {
                    println!(
                        "Reference coordinates: {}-{} (length={})",
                        ref_start, ref_end, ref_length
                    );
                    println!("Expected from issue: 3367501-3367604 (length=103)");
                    println!(
                        "Current result: {}-{} (length={})",
                        ref_start, ref_end, ref_length
                    );
                } else {
                    println!("Failed to lift coordinates");
                }
            }
            Err(e) => println!("Error lifting coordinates: {}", e),
        }

        break; // Only process first record
    }

    Ok(())
}

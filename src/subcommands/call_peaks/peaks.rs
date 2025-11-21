use anyhow::Result;

use super::fdr::FdrEntry;
use super::pileup::chrom_names_and_lengths;
use crate::cli::CallPeaksOptions;

/// Call peaks using the FDR table
///
/// This will:
/// 1. Run pileup on real data with FDR annotation
/// 2. Find local maxima
/// 3. Call peaks above FDR threshold
/// 4. Merge overlapping peaks
/// 5. Output results to BED file
pub fn call_peaks_with_fdr(
    _opts: &mut CallPeaksOptions,
    _bam: &mut rust_htslib::bam::IndexedReader,
    header: &rust_htslib::bam::HeaderView,
    fdr_table: &[FdrEntry],
) -> Result<()> {
    log::info!(
        "Calling peaks using FDR table with {} entries",
        fdr_table.len()
    );

    log::warn!("Peak calling logic not yet implemented - only FDR table generation is complete");
    for (chrom_str, chrom_len) in chrom_names_and_lengths(header)? {
        log::info!(
            "Finding peaks on chromosome {} with length {}",
            chrom_str,
            chrom_len
        );

        // Process fibers to build the track
        // Calculate scores and rolling max
        // TODO: Find local maxima
        // TODO: Annotate with FDR
        // TODO: Call peaks
    }

    Ok(())
}

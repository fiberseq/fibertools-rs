mod fdr;
mod peaks;
mod pileup;

pub use fdr::{make_fdr_table, write_fdr_table, FdrEntry, PileupRecord};
pub use peaks::call_peaks_with_fdr;
pub use pileup::{chrom_names_and_lengths, fibers_from_chromosome, process_chromosome_pileup_both};

use crate::cli::CallPeaksOptions;
use anyhow::Result;

pub fn run_call_peaks(opts: &mut CallPeaksOptions) -> Result<()> {
    log::info!("Starting FIRE peak calling");
    log::info!("  Input BAM: {}", opts.input.bam);
    log::info!("  Output: {}", opts.out);
    log::info!("  Max FDR: {}", opts.max_fdr);
    log::info!("  Window size: {}", opts.window_size);

    let mut bam = opts.input.indexed_bam_reader();
    let header = opts.input.header_view();

    if let Some(ref shuffled) = opts.shuffled {
        log::info!("  Shuffled positions file: {}", shuffled);
    } else {
        log::info!("  Using random shuffling for FDR calculation");
    }

    // Generate or load FDR table
    let fdr_table = if let Some(ref fdr_table_path) = opts.fdr_table {
        log::info!("Loading FDR table from: {}", fdr_table_path);
        // TODO: load from file
        Vec::new()
    } else {
        // Generate pileup for both real and shuffled data
        log::info!("Running pileup for real and shuffled data...");
        let (real_pileup, shuffled_pileup) = pileup::generate_pileups(opts, &mut bam, &header)?;

        // Generate FDR table from pileup data
        make_fdr_table(real_pileup, shuffled_pileup, opts.max_fdr)?
    };

    // Write FDR table if requested
    if let Some(ref fdr_out) = opts.fdr_table_out {
        log::info!("Writing FDR table to: {}", fdr_out);
        write_fdr_table(&fdr_table, fdr_out)?;
    }

    // Call peaks using the FDR table
    call_peaks_with_fdr(opts, &mut bam, &header, &fdr_table)?;

    log::info!("FIRE peak calling completed");
    Ok(())
}

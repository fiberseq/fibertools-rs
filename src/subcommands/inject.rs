use crate::cli::InjectOptions;
use crate::utils::fibertig::FiberTig;
use anyhow::Result;

pub fn run_inject(opts: &InjectOptions) -> Result<()> {
    log::info!("Creating mock BAM from reference FASTA: {}", opts.reference);

    // Create FiberTig from FASTA
    let fibertig = FiberTig::from_inject_opts(opts)?;

    log::info!("Generated {} sequences", fibertig.records().len());

    // Write to BAM file (with built-in broken pipe handling)
    fibertig.write_to_bam(opts)?;

    log::info!("Mock BAM written to: {}", opts.out);

    Ok(())
}

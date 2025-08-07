use crate::cli::PgInjectOptions;
use crate::utils::fibertig::FiberTig;
use anyhow::Result;

pub fn run_pg_inject(opts: &PgInjectOptions) -> Result<()> {
    if opts.extract {
        // Extract mode: read from BAM and output BED
        log::info!("Extracting BED annotations from BAM: {}", opts.reference);
        FiberTig::extract_to_bed(opts)?;
        log::info!("BED annotations extracted to: {}", opts.out);
    } else {
        // Inject mode: read from FASTA and output BAM
        log::info!("Creating mock BAM from reference FASTA: {}", opts.reference);

        // Create FiberTig from FASTA
        let fibertig = FiberTig::from_inject_opts(opts)?;

        log::info!("Generated {} sequences", fibertig.records().len());

        // Write to BAM file (with built-in broken pipe handling)
        fibertig.write_to_bam(opts)?;

        log::info!("Mock BAM written to: {}", opts.out);
    }

    Ok(())
}

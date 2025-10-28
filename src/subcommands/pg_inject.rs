use crate::cli::PgInjectOptions;
use crate::subcommands::pg_pansn::{apply_pansn_transformations, haplotype};
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
        let mut fibertig = FiberTig::from_inject_opts(opts)?;

        log::info!("Generated {} sequences", fibertig.records().len());

        // Apply panSN transformations if any are specified
        if opts.pansn.has_operations() {
            apply_pansn_transformations(&mut fibertig.header, &opts.pansn)?;

            // Apply haplotype tagging to records if both haplotype tags are provided
            if opts.pansn.hap1_tag.is_some() && opts.pansn.hap2_tag.is_some() {
                let haplotype_map = haplotype::build_haplotype_map(
                    &fibertig.header,
                    opts.pansn.hap1_tag.as_deref(),
                    opts.pansn.hap2_tag.as_deref(),
                )?;

                if let Some(hap_map) = haplotype_map {
                    for record in &mut fibertig.records {
                        haplotype::add_haplotype_tag(record, &hap_map, opts.pansn.min_mapq)?;
                    }
                }
            }
        }

        // Write to BAM file (with built-in broken pipe handling)
        fibertig.write_to_bam(opts)?;

        log::info!("Mock BAM written to: {}", opts.out);
    }

    Ok(())
}

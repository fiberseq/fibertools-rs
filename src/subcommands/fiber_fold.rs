use crate::cli::{FiberFoldOptions, PileupOptions};
use crate::subcommands::pileup::FiberseqPileup;
use crate::*;
use rust_htslib::bam::{FetchDefinition, IndexedReader};

/// get a pileup using the indexed bam reader for a fetch region
/// Returns a FiberseqPileup instead of writing to file
pub fn fold_pileup<'a>(
    bam: &'a mut IndexedReader,
    pileup_opts: &'a PileupOptions,
    rgn: FetchDefinition,
) -> Result<FiberseqPileup<'a>> {
    let (chrom, chrom_start, mut chrom_end) = match &rgn {
        FetchDefinition::RegionString(chrom, start, end) => (chrom, *start, *end),
        FetchDefinition::String(chrom_bytes) => {
            let tid = bam.header().tid(chrom_bytes).unwrap();
            let chrom_len = bam.header().target_len(tid).unwrap() as i64;
            (chrom_bytes, 0, chrom_len)
        }
        _ => return Err(anyhow::anyhow!("Unsupported fetch definition")),
    };
    let tid = bam.header().tid(chrom).unwrap();
    let chrom_len = bam.header().target_len(tid).unwrap() as i64;
    // lossy chrom as str
    let chrom = String::from_utf8_lossy(chrom).to_string();

    if chrom_end > chrom_len {
        chrom_end = chrom_len;
    }

    // fetch the data
    bam.fetch(rgn)?;
    let records = bam.records();

    // make the pileup
    log::debug!(
        "Initializing pileup for {}:{}-{}",
        chrom,
        chrom_start,
        chrom_end
    );
    let mut pileup = FiberseqPileup::new(
        &chrom,
        chrom_start as usize,
        chrom_end as usize,
        pileup_opts,
        &None, // No shuffled fibers
    );
    pileup.add_records(records)?;

    Ok(pileup)
}

pub fn fiber_fold(opts: &mut FiberFoldOptions) -> Result<()> {
    // Create a new InputBam by cloning the necessary fields
    let input_bam = crate::utils::input_bam::InputBam {
        bam: opts.input.bam.clone(),
        filters: opts.input.filters.clone(),
        global: opts.input.global.clone(),
        header: None,
    };
    
    // define the options used for creating pileup(s)
    let mut pileup_opts = crate::cli::PileupOptions {
        input: input_bam,
        rgn: None,
        out: "-".to_string(),
        m6a: true,
        cpg: true,
        haps: true,
        keep_zeros: true,
        per_base: true,
        fiber_coverage: true,
        shuffle: None,
        rolling_max: None,
        no_msp: false,
        no_nuc: false,
    };
    
    // Now create the bam reader using the input from pileup_opts
    let mut bam = pileup_opts.input.indexed_bam_reader();


    // TODO in a loop over the genome call fold_pileup and process the results

    todo!("FiberFold functionality not yet implemented")
}

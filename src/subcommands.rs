/// Add nucleosomes to a bam file
pub mod add_nucleosomes;
/// Center fiberseq information around a reference position
pub mod center;
/// Clear HiFi kinetics tags from a bam file
pub mod clear_kinetics;
pub mod ddda_to_m6a;
/// add decorators
pub mod decorator;
/// Extract fiberseq data into plain text formats
pub mod extract;
/// add fire data
pub mod fire;
pub mod footprint;
/// make a fire track from a bam file
pub mod pileup;
/// m6A prediction
pub mod predict_m6a;
/// Collect QC metrics
pub mod qc;
/// Remove base modifications from a bam record
pub mod strip_basemods;

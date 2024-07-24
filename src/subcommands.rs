/// Center fiberseq information around a reference position
pub mod center;
pub mod ddda_to_m6a;
/// add decorators
pub mod decorator;
/// Extract fiberseq data into plain text formats
pub mod extract;
/// add fire data
pub mod fire;
pub mod footprint;
/// Add nucleosomes to a bam file
pub mod nucleosomes;
/// make a fire track from a bam file
pub mod pileup;
/// m6A prediction
pub mod predict_m6a;
/// Remove base modifications from a bam record
pub mod strip_basemods;

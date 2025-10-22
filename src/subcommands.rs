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

pub mod fiber_hmm;
/// Create mock BAM from reference FASTA
pub mod pg_inject;
/// Lift annotations through a pangenome graph from source to target coordinates
pub mod pg_lift;
/// Add or strip panSN-spec prefixes from BAM contig names
pub mod pg_pansn;

pub mod validate;

//
pub mod merge_peaks;

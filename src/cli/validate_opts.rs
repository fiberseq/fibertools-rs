use crate::utils::input_bam::InputBam;
use clap::Args;

#[derive(Args, Debug)]
pub struct ValidateOptions {
    /// Path to the input BAM file
    #[clap(flatten)]
    pub bam: InputBam,

    /// Number of reads to validate
    #[clap(short, long, default_value = "1000")]
    pub num_reads: usize,

    /// Check for FIRE calls in the reads
    #[clap(short = 'f', long)]
    pub check_fire: bool,

    /// The fraction of reads that must be valid for the BAM to be considered valid
    #[clap(short, long, default_value = "0.5")]
    pub min_valid_fraction: f64,
}

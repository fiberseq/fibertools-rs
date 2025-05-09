use crate::utils::input_bam::InputBam;
use clap::Args;

#[derive(Args, Debug)]
pub struct ValidateOptions {
    /// Path to the input BAM file
    #[clap(flatten)]
    pub input_file: InputBam,

    /// Number of reads to validate
    #[clap(
        short = 'n',
        long = "num-reads",
        default_value = "1000",
        help = "Number of reads to validate"
    )]
    pub num_reads: usize,

    /// Check for FIRE calls in the reads
    #[clap(
        short = 'f',
        long = "check-fire",
        help = "Check for FIRE calls in the reads"
    )]
    pub check_fire: bool,
}

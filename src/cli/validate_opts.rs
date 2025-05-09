use crate::utils::input_bam::InputBam;
use clap::Args;

#[derive(Args, Debug)]
pub struct ValidateOptions {
    #[clap(flatten)]
    pub bam: InputBam,

    /// Number of reads to validate
    #[clap(short, long, default_value = "5000")]
    pub reads: usize,

    /// The fraction of reads that must have m6A calls to pass validation
    #[clap(short, long, default_value = "0.5")]
    pub m6a: f64,

    /// The fraction of reads that must have nucleosome and MSP calls to pass validation
    #[clap(short = 'n', long, default_value = "0.5")]
    pub nuc: f64,

    /// Check for FIRE calls in the reads, there must be at least one FIRE call to pass validation.
    #[clap(short = 'f', long)]
    pub fire: bool,

    /// Check for the fraction of reads with alignment to a reference genome
    #[clap(short, long, default_value = "0.0")]
    pub aligned: f64,

    /// Check for the fraction of reads with phasing information
    #[clap(short, long, default_value = "0.0")]
    pub phased: f64,

    /// Check for the fraction of reads with kinetics information
    #[clap(short, long, default_value = "0.0")]
    pub kinetics: f64,
}

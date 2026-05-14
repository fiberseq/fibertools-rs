use super::nucleosome_opts::NucleosomeParameters;
use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct PredictM6AOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output bam file with m6A calls in new/extended MM and ML bam tags
    #[clap(default_value = "-")]
    pub out: String,
    #[clap(flatten)]
    pub nuc: NucleosomeParameters,
    /// Keep hifi kinetics data
    #[clap(short, long)]
    pub keep: bool,
    /// Force a different minimum ML score
    #[clap(long, help_heading = "Developer-Options")]
    pub force_min_ml_score: Option<u8>,
    /// Keep all m6A calls regardless of how low the ML value is
    #[clap(long, help_heading = "Developer-Options")]
    pub all_calls: bool,
    /// Number of reads to include in batch prediction
    ///
    /// Increasing improves GPU performance at the cost of memory.
    #[clap(short, long, default_value = "1", help_heading = "Developer-Options")]
    pub batch_size: usize,
    /// Skip the actual prediction step to allow for testing the speed of other parts of the code
    #[clap(long, help_heading = "Developer-Options", hide = true)]
    pub fake: bool,
}

impl std::default::Default for PredictM6AOptions {
    fn default() -> Self {
        Self {
            input: InputBam::default(),
            out: "-".to_string(),
            nuc: NucleosomeParameters::default(),
            keep: false,
            force_min_ml_score: None,
            all_calls: false,
            batch_size: 1,
            fake: false,
        }
    }
}

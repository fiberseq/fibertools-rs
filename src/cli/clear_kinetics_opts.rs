use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct ClearKineticsOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output bam file without hifi kinetics
    #[clap(default_value = "-")]
    pub out: String,
}

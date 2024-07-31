use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct DddaToM6aOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output bam file
    #[clap(default_value = "-")]
    pub out: String,
}

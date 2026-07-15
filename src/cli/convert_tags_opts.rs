use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct ConvertTagsOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output bam file with MA spec tags
    #[clap(default_value = "-")]
    pub out: String,
}

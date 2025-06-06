use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct FiberFoldOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output path for the FiberFold results
    #[clap(short, long)]
    pub output: Option<String>,
}

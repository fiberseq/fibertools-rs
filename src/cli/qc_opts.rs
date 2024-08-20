use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct QcOpts {
    #[clap(flatten)]
    pub input: InputBam,
    #[clap(short, long, default_value = "-")]
    pub out: String,
}

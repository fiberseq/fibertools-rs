use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;
#[derive(Args, Debug, Clone)]

pub struct FiberHmmParameters {
    /// Minium nucleosome length
    #[clap(short, long, default_value = "90")]
    pub nucleosome_length: i64,
}

#[derive(Args, Debug)]
pub struct FiberHmmOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output bam file with nucleosome calls
    #[clap(default_value = "-")]
    pub out: String,
    #[clap(flatten)]
    pub nuc: FiberHmmParameters,
}

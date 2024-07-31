use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct DecoratorOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output path for bed12 file to be decorated
    #[clap(short, long)]
    pub bed12: String,
    /// Output path for decorator bed file
    #[clap(short, long, default_value = "-")]
    pub decorator: String,
}

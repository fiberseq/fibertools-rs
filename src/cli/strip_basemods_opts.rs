use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct StripBasemodsOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output bam file
    #[clap(default_value = "-")]
    pub out: String,
    #[clap(short, long, value_parser(["m6A","6mA", "5mC","CpG"]))]
    /// base modification to strip out of the bam file
    pub basemod: String,
    /// Drop forward strand of base modifications
    #[clap(long)]
    pub drop_forward: bool,
    /// Drop reverse strand of base modifications
    #[clap(long)]
    pub drop_reverse: bool,
}

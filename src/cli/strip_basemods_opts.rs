use super::super::utils::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct StripBasemodsOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output bam file
    #[clap(default_value = "-")]
    pub out: String,
    #[clap(short, long, default_value = "m6A",  value_parser(["m6A","6mA", "5mC","CpG"]))]
    /// base modification to strip out of the bam file
    pub basemod: String,
}

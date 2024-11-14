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
    pub basemod: Option<String>,
    /// filter out m6A modifications with less than this ML value
    #[clap(long, default_value = "0")]
    pub ml_m6a: u8,
    /// filter out 5mC modifications with less than this ML value
    #[clap(long, default_value = "0")]
    pub ml_5mc: u8,
    /// Drop forward strand of base modifications
    #[clap(long)]
    pub drop_forward: bool,
    /// Drop reverse strand of base modifications
    #[clap(long)]
    pub drop_reverse: bool,
}

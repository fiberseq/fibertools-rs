use super::super::utils::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct PileupOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Region string to make a pileup of. e.g. chr1:1-1000 or chr1:1-1,000
    /// If not provided will make a pileup of the whole genome
    #[clap(default_value = None)]
    pub rgn: Option<String>,
    /// Output file
    #[clap(short, long, default_value = "-")]
    pub out: String,
    /// include m6A calls
    #[clap(short, long)]
    pub m6a: bool,
    /// include 5mC calls
    #[clap(short, long)]
    pub cpg: bool,
    /// For each column add two new columns with the hap1 and hap2 specific data.
    #[clap(long)]
    pub haps: bool,
    /// Keep zero coverage regions
    #[clap(short, long)]
    pub keep_zeros: bool,
    /// Write output one base at a time even if the values do not change
    #[clap(short, long)]
    pub per_base: bool,
}

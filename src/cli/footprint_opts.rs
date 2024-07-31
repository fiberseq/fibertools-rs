use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct FootprintOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// BED file with the regions to footprint. Should all contain the same motif with proper strand information, and ideally be ChIP-seq peaks.
    #[clap(short, long)]
    pub bed: String,
    /// yaml describing the modules of the footprint
    #[clap(short, long)]
    pub yaml: String,
    /// Output bam
    #[clap(short, long, default_value = "-")]
    pub out: String,
}

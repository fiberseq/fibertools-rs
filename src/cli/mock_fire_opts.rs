use crate::cli::GlobalOpts;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct MockFireOptions {
    /// Input BED file where intervals become FIRE elements.
    /// The 4th column (name) groups intervals into the same mock read.
    #[clap()]
    pub bed: String,
    /// Output BAM file
    #[clap(short, long, default_value = "-")]
    pub out: String,
    /// Length of mock reads (default: auto-calculated from BED intervals)
    #[clap(short, long)]
    pub read_length: Option<i64>,
    /// FIRE quality score to assign (0-255, higher = more confident FIRE call)
    #[clap(short, long, default_value_t = 255)]
    pub quality: u8,
    /// Uncompressed BAM output (default: compressed)
    #[clap(short, long)]
    pub uncompressed: bool,
    #[clap(flatten)]
    pub global: GlobalOpts,
}

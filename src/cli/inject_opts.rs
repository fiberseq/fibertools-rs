use crate::cli::GlobalOpts;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct InjectOptions {
    /// Reference FASTA file to create mock BAM from (supports .gz/.bgz compression)
    #[clap()]
    pub reference: String,
    /// Output BAM file
    #[clap(short, long, default_value = "-")]
    pub out: String,
    /// Split contigs into multiple BAM records every N base pairs (0 = no splitting)
    #[clap(long, default_value = "10_000_000")]
    pub split_size: usize,
    #[clap(flatten)]
    pub global: GlobalOpts,
}

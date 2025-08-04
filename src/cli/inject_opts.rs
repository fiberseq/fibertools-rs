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
    /// Prefix to add to the contig names in the BAM header (PanSN-spec format)
    #[clap(long, short)]
    pub pansn_prefix: Option<String>,
    /// Split contigs into multiple BAM records every N base pairs (0 = no splitting)
    #[clap(short, long, default_value_t = 10_000_000)]
    pub split_size: usize,
    #[clap(flatten)]
    pub global: GlobalOpts,
}

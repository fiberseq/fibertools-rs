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
}

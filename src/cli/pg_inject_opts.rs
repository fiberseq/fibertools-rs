use crate::cli::{GlobalOpts, PansnParameters};
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct PgInjectOptions {
    /// Reference FASTA file to create mock BAM from (supports .gz/.bgz compression)
    #[clap()]
    pub reference: String,
    /// Output BAM file
    #[clap(short, long, default_value = "-")]
    pub out: String,
    /// Split contigs into multiple BAM records every N base pairs (0 = no splitting)
    #[clap(short, long, default_value_t = 50_000)]
    pub split_size: usize,
    /// Uncompressed BAM output (default: compressed)
    #[clap(short, long)]
    pub uncompressed: bool,
    /// Optional BED file with annotations to add to mock BAM records
    #[clap(short, long)]
    pub bed: Option<String>,
    /// Extract BED annotations from an annotated BAM file (reverses injection)
    #[clap(short, long)]
    pub extract: bool,
    /// Additionally write the output BAM header to this file
    #[clap(long)]
    pub header_out: Option<String>,
    #[clap(flatten)]
    pub global: GlobalOpts,
    #[clap(flatten)]
    pub pansn: PansnParameters,
}

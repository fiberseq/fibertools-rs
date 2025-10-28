use crate::cli::GlobalOpts;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct PgLiftOptions {
    /// Reference FASTA file to use as source coordinates
    #[clap()]
    pub reference: String,
    /// Input file with annotations to lift (defaults to BED format)
    #[clap(short, long)]
    pub input: String,
    /// Input file is in BAM format (default: BED format)
    #[clap(long)]
    pub bam: bool,
    /// Pangenome graph file in GBZ format
    #[clap(short, long)]
    pub graph: String,
    /// Output BED file with lifted annotations
    #[clap(short, long, default_value = "-")]
    pub out: String,
    /// panSN-spec prefix to add before graph injection (e.g., "HG002_2#0#")
    #[clap(short, long)]
    pub prefix: String,
    /// Target sample/haplotype name for surjection (e.g., "HG002_1")
    #[clap(long)]
    pub target: String,
    /// Split contigs into multiple BAM records every N base pairs for graph injection
    #[clap(long, default_value_t = 100_000)]
    pub split_size: usize,
    /// Number of threads to use for vg operations
    #[clap(long, default_value_t = 16)]
    pub vg_threads: usize,
    /// Path to vg binary (default: searches PATH)
    #[clap(long, default_value = "vg")]
    pub vg_binary: String,
    /// Delimiter character to use when stripping panSN-spec (default: '#')
    #[clap(short, long, default_value = "#")]
    pub delimiter: char,
    /// Keep intermediate files for debugging
    #[clap(long)]
    pub keep_intermediate: bool,
    #[clap(flatten)]
    pub global: GlobalOpts,
}

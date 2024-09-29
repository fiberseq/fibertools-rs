use crate::utils::input_bam::InputBam;
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
    /// Calculate coverage starting from the first MSP/NUC to the last MSP/NUC
    /// position instead of the complete span of the read alignment.
    #[clap(long)]
    pub fiber_coverage: bool,
    /// Shuffle the fiber-seq data according to a bed file of
    /// the shuffled positions of the fiber-seq data
    ///
    /// The bed file should have the following format:
    /// #chrom shuffled_start shuffled_end read_name original_start
    #[clap(long)]
    pub shuffle: Option<String>,
    /// Output a rolling max of the score column over X bases
    #[clap(long)]
    pub rolling_max: Option<usize>,
    /// No MSP columns
    #[clap(long)]
    pub no_msp: bool,
    /// No NUC columns
    #[clap(long)]
    pub no_nuc: bool,
}

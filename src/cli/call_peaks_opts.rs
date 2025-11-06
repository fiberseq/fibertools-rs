use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct CallPeaksOptions {
    #[clap(flatten)]
    pub input: InputBam,

    /// BED file with shuffled fiber positions (from bedtools shuffle)
    /// If not provided, will use all positions as real data (no FDR calculation)
    #[clap(short, long)]
    pub shuffled: Option<String>,

    /// Output BED file with called peaks
    #[clap(short, long, default_value = "-")]
    pub out: String,

    /// Maximum coverage threshold for filtering (optional)
    #[clap(long)]
    pub max_cov: Option<u32>,

    /// Minimum coverage threshold for filtering (optional)
    #[clap(long)]
    pub min_cov: Option<u32>,

    /// Maximum FDR threshold for peak calling
    #[clap(long, default_value = "0.05")]
    pub max_fdr: f64,

    /// Minimum fraction of accessible bases in peak
    #[clap(long, default_value = "0.0")]
    pub min_frac_accessible: f64,

    /// Rolling window size for finding local maxima (in base pairs)
    #[clap(long, default_value = "200")]
    pub window_size: usize,

    /// Minimum fraction of overlapping FIRE elements for merging peaks
    #[clap(long, default_value = "0.5")]
    pub min_frac_overlap: f64,

    /// Minimum reciprocal overlap for merging peaks
    #[clap(long, default_value = "0.75")]
    pub min_reciprocal_overlap: f64,

    /// Maximum number of grouping iterations for merging
    #[clap(long, default_value = "10")]
    pub max_grouping_iterations: usize,

    /// Skip the FDR table generation and use existing table
    #[clap(long)]
    pub fdr_table: Option<String>,

    /// Output the FDR table to this file
    #[clap(long)]
    pub fdr_table_out: Option<String>,

    /// Include nucleosome and MSP coverage in pileup (default: only FIRE coverage)
    #[clap(long)]
    pub include_nuc_msp: bool,

    /// Include haplotype-specific calls
    #[clap(long)]
    pub haps: bool,
}

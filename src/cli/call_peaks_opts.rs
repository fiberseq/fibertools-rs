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
    pub max_cov: Option<i32>,

    /// Minimum coverage threshold for filtering (optional)
    #[clap(long)]
    pub min_cov: Option<i32>,

    /// Number of standard deviations from median coverage to use for filtering (default: 5)
    /// If set, will calculate median +/- (sd_cov * std_dev) and use those as min/max coverage
    /// This overrides --max-cov and --min-cov if those are not explicitly set
    #[clap(long, default_value = "5.0")]
    pub sd_cov: f64,

    /// Maximum FDR threshold for peak calling (ignored if --min-fire-frac is set)
    #[clap(long, default_value = "0.05")]
    pub max_fdr: f64,

    /// Minimum fraction of fibers with FIREs required to call a peak
    /// If set, skips FDR calculation and uses this threshold instead
    /// For example, 0.5 means at least 50% of fibers must have a FIRE at the peak position
    #[clap(long)]
    pub min_fire_frac: Option<f64>,

    /// Minimum fraction of fibers with FIREs required as an additional filter (applied WITH FDR)
    /// Unlike --min-fire-frac, this is applied in addition to FDR filtering, not instead of it
    /// For example, 0.3 means at least 30% of fibers must have a FIRE AND FDR must be <= max_fdr
    #[clap(long, default_value = "0.1")]
    pub min_fire_frac_filter: f64,

    /// Minimum fraction of accessible bases in peak
    #[clap(long, default_value = "0.0", hide = true)]
    pub min_frac_accessible: f64,

    /// Rolling window size for finding local maxima (in base pairs)
    #[clap(long, default_value = "200")]
    pub window_size: usize,

    /// Minimum fraction of overlapping FIRE elements for merging peaks (Phase 2)
    #[clap(long, default_value = "0.5")]
    pub min_frac_overlap: f64,

    /// Minimum reciprocal overlap for merging peaks (Phase 3)
    #[clap(long, default_value = "0.75")]
    pub min_reciprocal_overlap: f64,

    /// High reciprocal overlap threshold for initial merging (Phase 1)
    #[clap(long, default_value = "0.90")]
    pub high_reciprocal_overlap: f64,

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

    /// Minimum FIRE coverage required to calculate a score (default: 4)
    #[clap(long, default_value = "4", hide = true)]
    pub min_fire_coverage: i32,
}

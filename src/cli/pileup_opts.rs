use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct PileupOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Region string(s) to make a pileup of. e.g. chr1:1-1000 or chr1:1-1,000
    /// Can be specified multiple times for multiple regions.
    /// If not provided will make a pileup of the whole genome
    #[clap(short, long)]
    pub rgn: Vec<String>,
    /// BED file with regions to query. If the BED file has a name column (4th column),
    /// the name will be added to each output line for that region.
    #[clap(short, long, conflicts_with = "rgn")]
    pub bed: Option<String>,
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
    /// Convenience: apply the FIRE peak-calling pipeline's fiber-level filters
    /// (`--skip-no-m6a`, `--min-msp 10`, `--min-ave-msp-size 10`) and enable
    /// `--fiber-coverage`. Individual filter flags still override when both
    /// are set.
    #[clap(long)]
    pub fire_filter: bool,
    /// Drop fibers with no m6A calls (matches `ft fire` behavior). Off by default;
    /// `--fire-filter` turns this on unless explicitly set to `false`.
    #[clap(long)]
    pub skip_no_m6a: Option<bool>,
    /// Drop fibers with fewer than `N` MSP calls (matches `ft fire`).
    /// Off (0) by default; `--fire-filter` sets this to 10 unless overridden.
    #[clap(long)]
    pub min_msp: Option<usize>,
    /// Drop fibers whose average MSP size is below `N` (matches `ft fire`).
    /// Off (0) by default; `--fire-filter` sets this to 10 unless overridden.
    #[clap(long)]
    pub min_ave_msp_size: Option<i64>,
}

impl PileupOptions {
    pub fn fire_filters(&self) -> crate::utils::fire::FireFiberFilters {
        crate::utils::fire::FireFiberFilters {
            skip_no_m6a: self.skip_no_m6a.unwrap_or(self.fire_filter),
            min_msp: self
                .min_msp
                .unwrap_or(if self.fire_filter { 10 } else { 0 }),
            min_ave_msp_size: self
                .min_ave_msp_size
                .unwrap_or(if self.fire_filter { 10 } else { 0 }),
        }
    }

    /// `--fire-filter` bundles `--fiber-coverage` in addition to the three
    /// filters. Callers should use this in place of reading `fiber_coverage`
    /// directly so the bundle stays coherent.
    pub fn effective_fiber_coverage(&self) -> bool {
        self.fiber_coverage || self.fire_filter
    }
}

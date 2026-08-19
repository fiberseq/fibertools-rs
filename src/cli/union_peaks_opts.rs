use crate::cli::GlobalOpts;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct UnionPeaksOptions {
    /// Input BED files, one per sample.
    /// Every interval becomes a FIRE element on a mock fiber for that sample, and
    /// peaks are called across all the samples at once. Overlapping intervals within
    /// one file are merged first, so a file can add at most 1 to a peak's support.
    #[clap(required = true, num_args = 1..)]
    pub beds: Vec<String>,
    /// Output BED file with union peaks
    #[clap(short, long, default_value = "-")]
    pub out: String,
    /// Sample names, comma separated, one per input BED [default: input file basenames]
    #[clap(long, value_delimiter = ',')]
    pub names: Vec<String>,
    /// Minimum number of input BEDs that must overlap a peak for it to be reported
    #[clap(short = 'n', long, default_value_t = 1)]
    pub min_support: usize,
    /// Rolling window size for finding local maxima (in base pairs)
    #[clap(long, default_value_t = 200)]
    pub window_size: usize,
    #[clap(flatten)]
    pub global: GlobalOpts,
}

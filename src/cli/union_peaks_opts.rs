use crate::cli::GlobalOpts;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct UnionPeaksOptions {
    /// Input BED files, one per sample.
    /// Every interval becomes a FIRE element on a mock fiber for that sample, and
    /// peaks are called across all the samples at once. Overlapping or book-ended intervals within
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
    /// Minimum fraction of input BEDs that must overlap a peak for it to be
    /// reported (0-1). Applied together with --min-support.
    #[clap(long, value_parser = frac_in_range)]
    pub min_frac_support: Option<f64>,
    /// Rolling window size for finding local maxima (in base pairs).
    /// Only local maxima are kept, so at most one peak is reported per window.
    #[clap(long, default_value_t = 200)]
    pub window_size: usize,
    #[clap(flatten)]
    pub global: GlobalOpts,
}

fn frac_in_range(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|e| format!("{e}"))?;
    if (0.0..=1.0).contains(&v) {
        Ok(v)
    } else {
        Err("must be between 0 and 1".to_string())
    }
}

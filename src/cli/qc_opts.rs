use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct QcOpts {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output text file with QC metrics. The format is a tab-separated file with the following columns: "statistic\tvalue\tcount" where "statistic" is the name of the metric, "value" is the value of the metric, and "count" is the number of times the metric was observed.
    #[clap(default_value = "-")]
    pub out: String,
    /// Calculate the auto-correlation function of the m6A marks in the fiber-seq data.
    #[clap(long)]
    pub acf: bool,
    /// maximum lag for the ACF calculation
    #[clap(long, default_value = "250")]
    pub acf_max_lag: usize,
    /// Minimum number of m6A marks to use a read in the ACF calculation
    #[clap(long, default_value = "100")]
    pub acf_min_m6a: usize,
    /// maximum number of reads to use in the ACF calculation
    #[clap(long, default_value = "10000")]
    pub acf_max_reads: usize,
    /// After sampling the first "acf-max-reads" randomly sample one of every "acf-sample-rate" reads and replace one of the previous reads at random.
    #[clap(long, default_value = "100")]
    pub acf_sample_rate: f32,
    /// In the output include a measure of the number of m6A events per MSPs of a given size.
    /// The output format is: "m6a_per_msp_size\t{m6A count},{MSP size},{is a FIRE}\t{count}"
    /// e.g. "m6a_per_msp_size\t35,100,false\t100"
    #[clap(short, long)]
    pub m6a_per_msp: bool,
    /// Only process the first "n" reads in the input bam file.
    #[clap(long)]
    pub n_reads: Option<usize>,
}

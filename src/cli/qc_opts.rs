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
    /// In the output include a measure of the number of m6A events per MSPs of a given size.
    /// The output format is: "m6a_per_msp_size\t{m6A count},{MSP size},{is a FIRE}\t{count}"
    /// e.g. "m6a_per_msp_size\t35,100,false\t100"
    #[clap(short, long)]
    pub m6a_per_msp: bool,
}

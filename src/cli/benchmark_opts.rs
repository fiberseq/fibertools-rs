use crate::utils::input_bam::InputBam;
use clap::Args;

#[derive(Args, Debug)]
pub struct BenchmarkOptions {
    #[clap(flatten)]
    pub input: InputBam,

    /// Region to benchmark (format: chr:start-end)
    #[clap(short, long)]
    pub region: Option<String>,
}

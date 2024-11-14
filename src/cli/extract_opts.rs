use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct ExtractOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Report positions in reference sequence coordinates
    #[clap(short, long, default_value = "true",
          default_value_ifs([
              ("molecular", "true", "false"),
              ("molecular", "false", "true"),
          ]))
      ]
    pub reference: bool,
    /// Report positions in the molecular sequence coordinates
    #[clap(long, default_value = "false")]
    pub molecular: bool,
    /// Output path for m6a bed12
    #[clap(long)]
    pub m6a: Option<String>,
    /// Output path for 5mC (CpG, primrose) bed12
    #[clap(short, long)]
    pub cpg: Option<String>,
    /// Output path for methylation sensitive patch (msp) bed12
    #[clap(long)]
    pub msp: Option<String>,
    /// Output path for nucleosome bed12
    #[clap(short, long)]
    pub nuc: Option<String>,
    /// Output path for a tabular format including "all" fiberseq information in the bam
    #[clap(short, long)]
    pub all: Option<String>,
    /// Include per base quality scores in "fiber_qual"
    #[clap(short, long, help_heading = "All-Format-Options")]
    pub quality: bool,
    /// Simplify output by removing fiber sequence
    #[clap(short, long, help_heading = "All-Format-Options")]
    pub simplify: bool,
}

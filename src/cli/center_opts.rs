use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct CenterOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Bed file on which to center fiberseq reads. Data is adjusted to the start position of the bed file and corrected for strand if the strand is indicated in the 6th column of the bed file. The 4th column will also be checked for the strand but only after the 6th is.        
    /// If you include strand information in the 4th (or 6th) column it will orient data accordingly and use the end position of bed record instead of the start if on the minus strand. This means that profiles of motifs in both the forward and minus orientation will align to the same central position.
    #[clap(short, long)]
    pub bed: String,
    /// Set a maximum distance from the start of the motif to keep a feature
    #[clap(short, long)]
    pub dist: Option<i64>,
    /// Provide data in wide format, one row per read
    #[clap(short, long)]
    pub wide: bool,
    /// Return relative reference position instead of relative molecular position
    #[clap(short, long)]
    pub reference: bool,
    /// Replace the sequence output column with just "N".
    #[clap(short, long)]
    pub simplify: bool,
}

use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

pub static NUC_LEN: &str = "75";
pub static COMBO_NUC_LEN: &str = "100";
pub static MIN_DIST_ADDED: &str = "25";
pub static DIST_FROM_END: &str = "45";
pub static ALLOWED_SKIPS: &str = "-1";

#[derive(Args, Debug, Clone)]
pub struct NucleosomeParameters {
    /// Minium nucleosome length
    #[clap(short, long, default_value = NUC_LEN)]
    pub nucleosome_length: i64,
    /// Minium nucleosome length when combining over a single m6A
    #[clap(short, long, default_value = COMBO_NUC_LEN)]
    pub combined_nucleosome_length: i64,
    /// Minium distance needed to add to an already existing nuc by crossing an m6a
    #[clap(long, default_value = MIN_DIST_ADDED)]
    pub min_distance_added: i64,
    /// Minimum distance from the end of a fiber to call a nucleosome or MSP
    #[clap(short, long, default_value = DIST_FROM_END)]
    pub distance_from_end: i64,
    /// Most m6A events we can skip over to get to the nucleosome length when using D-segment algorithm. 2 is often a good value, negative values disable D-segment for the simple caller.
    #[clap(short, long, default_value = ALLOWED_SKIPS, hide = true)]
    pub allowed_m6a_skips: i64,
}

impl std::default::Default for NucleosomeParameters {
    fn default() -> Self {
        Self {
            nucleosome_length: NUC_LEN.parse().unwrap(),
            combined_nucleosome_length: COMBO_NUC_LEN.parse().unwrap(),
            min_distance_added: MIN_DIST_ADDED.parse().unwrap(),
            distance_from_end: DIST_FROM_END.parse().unwrap(),
            allowed_m6a_skips: ALLOWED_SKIPS.parse().unwrap(),
        }
    }
}

#[derive(Args, Debug)]
pub struct AddNucleosomeOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output bam file with nucleosome calls
    #[clap(default_value = "-")]
    pub out: String,
    #[clap(flatten)]
    pub nuc: NucleosomeParameters,
}

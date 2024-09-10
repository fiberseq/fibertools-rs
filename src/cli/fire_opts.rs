use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct FireOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output file (BAM by default, table of MSP features if `--feats-to-text` is used, and bed9 + if `--extract`` is used)
    #[clap(default_value = "-")]
    pub out: String,
    /// Use a ONT heuristic adjustment for FIRE calling.
    /// This adjusts the observed number of m6A counts by adding pseudo counts to account for the single stranded nature of ONT data.
    #[clap(long, env)]
    pub ont: bool,
    /// Use a human model for FIRE calling
    #[clap(hide = true, long, env)]
    pub human: bool,
    /// Use a yeast model for FIRE calling
    #[clap(hide = true, long, env)]
    pub yeast: bool,
    /// Output just FIRE elements in bed9 format
    #[clap(short, long)]
    pub extract: bool,
    /// When extracting bed9 format include all MSPs and nucleosomes
    #[clap(long)]
    pub all: bool,
    /// Output FIREs features for training in a table format
    #[clap(short, long)]
    pub feats_to_text: bool,
    /// Don't write reads with no m6A calls to the output bam
    #[clap(short, long)]
    pub skip_no_m6a: bool,
    /// Skip reads without at least `N` MSP calls
    #[clap(long, default_value = "0", env = "MIN_MSP")]
    pub min_msp: usize,
    /// Skip reads without an average MSP size greater than `N`
    #[clap(long, default_value = "0", env)]
    pub min_ave_msp_size: i64,
    /// Width of bin for feature collection
    #[clap(short, long, default_value = "40", env, 
        default_value_ifs([
            ("human", "true", "40"),
            ("yeast", "true", "100")
        ])
    )]
    pub width_bin: i64,
    /// Number of bins to collect
    #[clap(short, long, default_value = "9", env)]
    pub bin_num: i64,
    /// Calculate stats for the highest X bp window within each MSP
    /// Should be a fair amount higher than the expected linker length.
    #[clap(long, default_value = "100", env)]
    pub best_window_size: i64,
    /// Use 5mC data in FIREs
    #[clap(short, long, hide = true)]
    pub use_5mc: bool,
    /// Minium length of msp to call a FIRE
    #[clap(long, default_value = "85", env)]
    pub min_msp_length_for_positive_fire_call: i64,
    /// Optional path to a model json file.
    /// If not provided ft will use the default model (recommended).
    #[clap(long, env = "FIRE_MODEL")]
    pub model: Option<String>,
    /// Optional path to a FDR table
    #[clap(long, env)]
    pub fdr_table: Option<String>,
}

impl Default for FireOptions {
    fn default() -> Self {
        Self {
            input: Default::default(),
            out: "-".to_string(),
            ont: false,
            human: false,
            yeast: false,
            extract: false,
            all: false,
            feats_to_text: false,
            skip_no_m6a: false,
            min_msp: 0,
            min_ave_msp_size: 0,
            width_bin: 40,
            bin_num: 9,
            best_window_size: 100,
            use_5mc: false,
            min_msp_length_for_positive_fire_call: 85,
            model: None,
            fdr_table: None,
        }
    }
}

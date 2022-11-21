use clap::IntoApp;
use clap::{AppSettings, Parser, Subcommand};

#[derive(Parser, Debug)]
#[clap(
    author,
    version,
    about,
    propagate_version = true,
    subcommand_required = true,
    infer_subcommands = true,
    arg_required_else_help = true,
    help_expected = true
)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
pub struct Cli {
    /// Threads for decompression.
    #[clap(
        global = true,
        short,
        long,
        default_value_t = 8,
        help_heading = "GLOBAL"
    )]
    pub threads: usize,

    /// Logging level [-v: Info, -vv: Debug, -vvv: Trace].
    #[clap(
        global = true,
        short,
        long,
        parse(from_occurrences),
        help_heading = "DEBUG"
    )]
    pub verbose: usize,
    /// Turn of all logging.
    #[clap(global = true, long, help_heading = "DEBUG")]
    pub quiet: bool,

    #[clap(subcommand)]
    pub command: Option<Commands>,
}

///
/// This structure contains all the subcommands for fiberseq-rs and their help descriptions.
///
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Extract fiberseq data into plain text files.
    #[clap(visible_aliases = &["ex", "e"])]
    Extract {
        /// fiberseq bam file
        #[clap(default_value = "-")]
        bam: String,
        /// report in reference sequence coordinates.
        #[clap(short, long)]
        reference: bool,
        /// simplify output by remove fiber sequence
        #[clap(short, long)]
        simplify: bool,
        /// include per base quality scores in "fiber_qual"
        #[clap(short, long)]
        quality: bool,
        /// Minium score in the ML tag to include in the output.
        #[clap(short, long, default_value = "150")]
        min_ml_score: u8,
        /// Output path for m6a bed12.
        #[clap(long)]
        m6a: Option<String>,
        /// Output path for 5mC (CpG, primrose) bed12.
        #[clap(short, long)]
        cpg: Option<String>,
        /// Output path for methylation sensitive patch (msp) bed12.
        #[clap(long)]
        msp: Option<String>,
        /// Output path for nucleosome bed12.
        #[clap(short, long)]
        nuc: Option<String>,
        /// Output path for .
        #[clap(short, long)]
        all: Option<String>,
        /// Add the full floating point predictions of the ML model
        #[clap(short, long)]
        full_float: bool,
    },
    /// This command centers fiberseq data around given reference positions. This is useful for making aggregate m6a and CpG observations, as well as visualization of SVs.
    #[clap(visible_aliases = &["c", "ct"])]
    Center {
        /// fiberseq bam file, must be aligned and have an index.
        bam: String,
        /// Bed file on which to center fiberseq reads. Data is adjusted to the start position of the bed file and corrected for strand if a 4th strand column is included.
        ///
        /// If you include strand information in the 4th column it will orient data accordingly and use the end position of bed record instead of the start if on the minus strand. This means that profiles of motifs in both the forward and minus orientation will align to the same central position.
        bed: String,
        /// Minium score in the ML tag to include in the output.
        #[clap(short, long, default_value = "150")]
        min_ml_score: u8,
        /// Provide data in wide format, one row per read.
        #[clap(short, long)]
        wide: bool,
    },
    /// Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags.
    #[clap(visible_aliases = &["m6A", "m6a"])]
    PredictM6A {
        /// Bam HiFi file with kinetics.
        #[clap(default_value = "-")]
        bam: String,
        /// Output bam file with m6A calls in new/extended MM and ML bam tags.
        #[clap(default_value = "-")]
        out: String,
        /// keep hifi kinetics data
        #[clap(short, long)]
        keep: bool,
        /// use CNN model for prediction instead of XGB.
        #[clap(short, long)]
        cnn: bool,
        /// Add a bam tag (mp) with the full floating point predictions of the ML model
        #[clap(short, long)]
        full_float: bool,
    },
}

pub fn make_cli_parse() -> Cli {
    Cli::parse()
}

pub fn make_cli_app() -> clap::Command<'static> {
    Cli::command()
}

use clap::{Command, CommandFactory, Parser, Subcommand, ValueHint};
use clap_complete::{generate, Generator, Shell};
use std::io;

#[derive(Parser, Debug)]
#[clap(
    author,
    version,
    about,
    propagate_version = true,
    subcommand_required = true,
    infer_subcommands = true,
    arg_required_else_help = true
)]
#[command(version = super::LONG_VERSION)]
pub struct Cli {
    /// Threads
    #[clap(
        global = true,
        short,
        long,
        default_value_t = 8,
        help_heading = "Global-Options"
    )]
    pub threads: usize,
    /// Logging level [-v: Info, -vv: Debug, -vvv: Trace]
    #[clap(
        global = true,
        short,
        long,
        action = clap::ArgAction::Count,
        help_heading = "Debug-Options"
    )]
    pub verbose: u8,
    /// Turn of all logging
    #[clap(global = true, long, help_heading = "Debug-Options")]
    pub quiet: bool,
    /// Subcommands for fibertools-rs
    #[clap(subcommand)]
    pub command: Option<Commands>,
}

///
/// This structure contains all the subcommands for fiberseq-rs and their help descriptions.
///
#[derive(Subcommand, Debug, PartialEq, Eq)]
pub enum Commands {
    /// Extract fiberseq data into plain text files
    #[clap(visible_aliases = &["ex", "e"])]
    Extract {
        /// Fiberseq bam file
        #[arg(default_value = "-", value_hint = ValueHint::AnyPath)]
        bam: String,
        /// Report in reference sequence coordinates
        #[clap(short, long)]
        reference: bool,
        /// Minium score in the ML tag to include in the output
        #[clap(short, long, default_value = "150")]
        min_ml_score: u8,
        /// Output path for m6a bed12
        #[clap(long)]
        m6a: Option<String>,
        /// Output path for 5mC (CpG, primrose) bed12
        #[clap(short, long)]
        cpg: Option<String>,
        /// Output path for methylation sensitive patch (msp) bed12
        #[clap(long)]
        msp: Option<String>,
        /// Output path for nucleosome bed12
        #[clap(short, long)]
        nuc: Option<String>,
        /// Output path for a tabular format including "all" fiberseq information in the bam
        #[clap(short, long)]
        all: Option<String>,
        /// Include per base quality scores in "fiber_qual"
        #[clap(short, long, help_heading = "All-Format-Options")]
        quality: bool,
        /// Add the full floating point predictions of the ML model
        #[clap(short, long, help_heading = "All-Format-Options")]
        full_float: bool,
        /// Simplify output by remove fiber sequence
        #[clap(short, long, help_heading = "All-Format-Options")]
        simplify: bool,
    },
    /// This command centers fiberseq data around given reference positions. This is useful for making aggregate m6A and CpG observations, as well as visualization of SVs
    #[clap(visible_aliases = &["c", "ct"])]
    Center {
        /// Fiberseq bam file, must be aligned and have an index
        bam: String,
        /// Bed file on which to center fiberseq reads. Data is adjusted to the start position of the bed file and corrected for strand if the strand is indicated in the 6th column of the bed file. The 4th column will also be checked for the strand but only after the 6th is.
        ///
        /// If you include strand information in the 4th column it will orient data accordingly and use the end position of bed record instead of the start if on the minus strand. This means that profiles of motifs in both the forward and minus orientation will align to the same central position.
        bed: String,
        /// Minium score in the ML tag to include in the output
        #[clap(short, long, default_value = "150")]
        min_ml_score: u8,
        /// Set a maximum distance from the start of the motif to keep a feature
        #[clap(short, long)]
        dist: Option<i64>,
        /// Provide data in wide format, one row per read
        #[clap(short, long)]
        wide: bool,
        /// Return relative reference position instead of relative molecular position
        #[clap(short, long)]
        reference: bool,
    },
    #[cfg(feature = "predict")]
    /// Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags
    #[clap(visible_aliases = &["m6A", "m6a"])]
    PredictM6A {
        /// Bam HiFi file with kinetics
        #[clap(default_value = "-")]
        bam: String,
        /// Output bam file with m6A calls in new/extended MM and ML bam tags
        #[clap(default_value = "-")]
        out: String,
        /// Minium nucleosome length
        #[clap(short, long, default_value = "75")]
        nucleosome_length: i64,
        /// Minium nucleosome length when combining over a single m6A
        #[clap(short, long, default_value = "100")]
        combined_nucleosome_length: i64,
        /// Minimum distance from the end of a fiber to call a nucleosome or MSP
        #[clap(short, long, default_value = "45")]
        distance_from_end: i64,
        /// Keep hifi kinetics data
        #[clap(short, long)]
        keep: bool,
        /// Set a minimum ML score to keep on instead of using the model specific minimum ML score.
        #[clap(short, long, help_heading = "Developer-Options")]
        min_ml_score: Option<u8>,
        /// Keep all m6A calls regardless of how low the ML value is
        #[clap(short, long, help_heading = "Developer-Options")]
        all_calls: bool,
        /// Use the XGBoost model for prediction
        #[clap(short, long)]
        xgb: bool,
        /// Use the CNN model for prediction
        #[clap(short, long, help_heading = "Developer-Options")]
        cnn: bool,
        /// Use the semi-supervised CNN model for prediction [default: true]
        #[clap(
            short,
            long,
            default_value_ifs([
                ("cnn", "true", "false"),
                ("xgb", "false", "true"),
            ]),
        )]
        semi: bool,
        /// Add a bam tag (mp) with the full floating point predictions of the ML model
        ///
        /// For debugging only.
        #[clap(short, long, help_heading = "Developer-Options")]
        full_float: bool,
        /// Number of reads to include in batch prediction
        ///
        /// Increasing improves GPU performance at the cost of memory.
        #[clap(
            short,
            long,
            default_value = "1",
            default_value_if("cnn", "true", "1"),
            help_heading = "Developer-Options"
        )]
        batch_size: usize,
    },
    /// Remove HiFi kinetics tags from the input bam file
    ClearKinetics {
        /// Bam HiFi file with kinetics
        #[clap(default_value = "-")]
        bam: String,
        /// Output bam file without hifi kinetics
        #[clap(default_value = "-")]
        out: String,
    },
    /// Add nucleosomes to a bam file with m6a predictions
    AddNucleosomes {
        /// Bam HiFi file with m6A calls
        #[clap(default_value = "-")]
        bam: String,
        /// Output bam file with nucleosome calls
        #[clap(default_value = "-")]
        out: String,
        /// Minium nucleosome length
        #[clap(short, long, default_value = "75")]
        nucleosome_length: i64,
        /// Minium nucleosome length when combining over a single m6A
        #[clap(short, long, default_value = "100")]
        combined_nucleosome_length: i64,
        /// Minimum distance from the end of a fiber to call a nucleosome or MSP
        #[clap(short, long, default_value = "45")]
        distance_from_end: i64,
    },
    /// Make command line completions
    Completions {
        /// If provided, outputs the completion file for given shell
        #[arg(value_enum)]
        shell: Shell,
    },
    /// Make a man page for fibertools-rs
    ///
    /// Writes file for `man` to stdout.
    Man {},
}

pub fn print_completions<G: Generator>(gen: G, cmd: &mut Command) {
    generate(gen, cmd, cmd.get_name().to_string(), &mut io::stdout());
}

pub fn make_cli_parse() -> Cli {
    Cli::parse()
}

pub fn make_cli_app() -> Command {
    Cli::command()
}

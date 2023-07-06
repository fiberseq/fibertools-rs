use super::nucleosomes::*;
use anstyle;
use clap::{Command, CommandFactory, Parser, Subcommand, ValueHint};
use clap_complete::{generate, Generator, Shell};
use std::io;

pub static MIN_ML_SCORE: &str = "125";

pub fn get_styles() -> clap::builder::Styles {
    let cmd_color = anstyle::AnsiColor::Magenta;
    let header_color = anstyle::AnsiColor::BrightGreen;
    let placeholder_color = anstyle::AnsiColor::Cyan;
    let header_style = anstyle::Style::new()
        .bold()
        .underline()
        .fg_color(Some(anstyle::Color::Ansi(header_color)));
    let cmd_style = anstyle::Style::new()
        .bold()
        .fg_color(Some(anstyle::Color::Ansi(cmd_color)));
    let placeholder_style = anstyle::Style::new()
        //.bold()
        .fg_color(Some(anstyle::Color::Ansi(placeholder_color)));
    clap::builder::Styles::styled()
        .header(header_style)
        .literal(cmd_style)
        .usage(header_style)
        .placeholder(placeholder_style)
}

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
#[command(styles=get_styles())]
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
    #[cfg(feature = "predict")]
    /// Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags. Also adds nucleosome (nl, ns) and MTase sensitive patches (al, as).
    #[clap(visible_aliases = &["m6A", "m6a"])]
    PredictM6A {
        /// Bam HiFi file with kinetics
        #[clap(default_value = "-")]
        bam: String,
        /// Output bam file with m6A calls in new/extended MM and ML bam tags
        #[clap(default_value = "-")]
        out: String,
        /// Minium nucleosome length
        #[clap(short, long, default_value = NUC_LEN)]
        nucleosome_length: i64,
        /// Minium nucleosome length when combining over a single m6A
        #[clap(short, long, default_value = COMBO_NUC_LEN)]
        combined_nucleosome_length: i64,
        /// Minium distance needed to add to an already existing nuc by crossing an m6a
        #[clap(long, default_value = MIN_DIST_ADDED)]
        min_distance_added: i64,
        /// Minimum distance from the end of a fiber to call a nucleosome or MSP
        #[clap(short, long, default_value = DIST_FROM_END)]
        distance_from_end: i64,
        /// Most m6A events we can skip over to get to the nucleosome length when using D-segment algorithm. 2 is often a good value, negative values disable D-segment for the simple caller.
        #[clap(long, default_value = ALLOWED_SKIPS, hide = true)]
        allowed_m6a_skips: i64,
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
        #[clap(long, help_heading = "Developer-Options")]
        xgb: bool,
        /// Use the CNN model for prediction
        #[clap(long, help_heading = "Developer-Options")]
        cnn: bool,
        /// Use the semi-supervised CNN model for prediction [default: true]
        #[clap(
            short,
            long,
            default_value_ifs([
                ("cnn", "true", "false"),
                ("xgb", "false", "true"),
            ]),
            help_heading = "Developer-Options",
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
    /// Add nucleosomes to a bam file with m6a predictions
    AddNucleosomes {
        /// Bam HiFi file with m6A calls
        #[clap(default_value = "-")]
        bam: String,
        /// Output bam file with nucleosome calls
        #[clap(default_value = "-")]
        out: String,
        /// Minium nucleosome length
        #[clap(short, long, default_value = NUC_LEN)]
        nucleosome_length: i64,
        /// Minium nucleosome length when combining over a single m6A
        #[clap(short, long, default_value = COMBO_NUC_LEN)]
        combined_nucleosome_length: i64,
        /// Minium distance needed to add to an already existing nuc by crossing an m6a
        #[clap(short, long, default_value = MIN_DIST_ADDED)]
        min_distance_added: i64,
        /// Minimum distance from the end of a fiber to call a nucleosome or MSP
        #[clap(short, long, default_value = DIST_FROM_END)]
        distance_from_end: i64,
        /// Most m6A events we can skip over to get to the nucleosome length when using D-segment algorithm. 2 is often a good value, negative values disable D-segment for the simple caller.
        #[clap(short, long, default_value = ALLOWED_SKIPS, hide = true)]
        allowed_m6a_skips: i64,
    },
    /// Extract fiberseq data into plain text files
    #[clap(visible_aliases = &["ex", "e"])]
    Extract {
        /// Fiberseq bam file
        #[arg(default_value = "-", value_hint = ValueHint::AnyPath)]
        bam: String,
        /// Report in reference sequence coordinates
        #[clap(short, long, default_value = "true",
            default_value_ifs([
                ("molecular", "true", "false"),
                ("molecular", "false", "true"),
            ]))
        ]
        reference: bool,
        /// Report positions in the molecular sequence coordinates
        #[clap(long, default_value = "false")]
        molecular: bool,
        /// Minium score in the ML tag to include in the output
        #[clap(short, long, default_value = MIN_ML_SCORE)]
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
        /// If you include strand information in the 4th (or 6th) column it will orient data accordingly and use the end position of bed record instead of the start if on the minus strand. This means that profiles of motifs in both the forward and minus orientation will align to the same central position.
        bed: String,
        /// Minium score in the ML tag to include in the output
        #[clap(short, long, default_value = MIN_ML_SCORE)]
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
    /// Remove HiFi kinetics tags from the input bam file
    ClearKinetics {
        /// Bam HiFi file with kinetics
        #[clap(default_value = "-")]
        bam: String,
        /// Output bam file without hifi kinetics
        #[clap(default_value = "-")]
        out: String,
    },
    /// Strip out select base modifications
    StripBasemods {
        /// Bam HiFi file with base mods
        #[clap(default_value = "-")]
        bam: String,
        /// Output bam file
        #[clap(default_value = "-")]
        out: String,
        #[clap(short, long, default_value = "m6A",  value_parser(["m6A","6mA", "5mC","CpG"]))]
        /// base modification to strip out of the bam file
        basemod: String,
    },
    /// Make command line completions
    #[clap(hide = true)]
    Completions {
        /// If provided, outputs the completion file for given shell
        #[arg(value_enum)]
        shell: Shell,
    },
    /// Make a man page for fibertools-rs
    ///
    /// Writes file for `man` to stdout.
    #[clap(hide = true)]
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

use super::nucleosomes::*;
use anstyle;
use clap::{Args, Command, CommandFactory, Parser, Subcommand, ValueHint};
use clap_complete::{generate, Generator, Shell};
use std::{fmt::Debug, io};

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

#[derive(Parser)]
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
    /// Turn off all logging
    #[clap(global = true, long, help_heading = "Debug-Options")]
    pub quiet: bool,
    /// Subcommands for fibertools-rs
    #[clap(subcommand)]
    pub command: Option<Commands>,
}

///
/// This structure contains all the subcommands for fiberseq-rs and their help descriptions.
///
#[derive(Subcommand)]
pub enum Commands {
    /// Predict m6A positions using HiFi kinetics data and encode the results in the MM and ML bam tags. Also adds nucleosome (nl, ns) and MTase sensitive patches (al, as).
    #[clap(visible_aliases = &["m6A", "m6a"])]
    PredictM6A(PredictM6AOptions),
    /// Add nucleosomes to a bam file with m6a predictions
    AddNucleosomes(AddNucleosomeOptions),
    /// Add FIREs (Fiber-seq Inferred Regulatory Elements) to a bam file with m6a predictions
    Fire(FireOptions),
    /// Extract fiberseq data into plain text files.
    ///
    /// See https://fiberseq.github.io/fibertools-rs/docs/extract.html for a description of the outputs.
    #[clap(visible_aliases = &["ex", "e"])]
    Extract(ExtractOptions),
    /// This command centers fiberseq data around given reference positions. This is useful for making aggregate m6A and CpG observations, as well as visualization of SVs.
    ///
    ///  See https://fiberseq.github.io/fibertools-rs/docs/center.html for a description of the output.
    #[clap(visible_aliases = &["c", "ct"])]
    Center(CenterOptions),
    /// Infer footprints from fiberseq data
    Footprint(FootprintOptions),
    /// Make decorated bed files for fiberseq data
    TrackDecorators(DecoratorOptions),
    /// make a pileup track from a FIRE bam
    Pileup(PileupOptions),
    /// Remove HiFi kinetics tags from the input bam file
    ClearKinetics(ClearKineticsOptions),
    /// Strip out select base modifications
    StripBasemods(StripBasemodsOptions),
    /// Make command line completions
    #[clap(hide = true)]
    Completions(CompletionOptions),
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

#[derive(Args, Debug, PartialEq, Eq)]
pub struct PredictM6AOptions {
    /// Bam HiFi file with kinetics
    #[clap(default_value = "-")]
    pub bam: String,
    /// Output bam file with m6A calls in new/extended MM and ML bam tags
    #[clap(default_value = "-")]
    pub out: String,
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
    #[clap(long, default_value = ALLOWED_SKIPS, hide = true)]
    pub allowed_m6a_skips: i64,
    /// Keep hifi kinetics data
    #[clap(short, long)]
    pub keep: bool,
    /// Set a minimum ML score to keep on instead of using the model specific minimum ML score.
    #[clap(short, long, help_heading = "Developer-Options")]
    pub min_ml_score: Option<u8>,
    /// Keep all m6A calls regardless of how low the ML value is
    #[clap(short, long, help_heading = "Developer-Options")]
    pub all_calls: bool,
    /// Number of reads to include in batch prediction
    ///
    /// Increasing improves GPU performance at the cost of memory.
    #[clap(short, long, default_value = "1", help_heading = "Developer-Options")]
    pub batch_size: usize,
}

impl std::default::Default for PredictM6AOptions {
    fn default() -> Self {
        Self {
            bam: "-".to_string(),
            out: "-".to_string(),
            nucleosome_length: NUC_LEN.parse().unwrap(),
            combined_nucleosome_length: COMBO_NUC_LEN.parse().unwrap(),
            min_distance_added: MIN_DIST_ADDED.parse().unwrap(),
            distance_from_end: DIST_FROM_END.parse().unwrap(),
            allowed_m6a_skips: ALLOWED_SKIPS.parse().unwrap(),
            keep: false,
            min_ml_score: None,
            all_calls: false,
            batch_size: 1,
        }
    }
}

#[derive(Args, Debug, PartialEq, Eq)]
pub struct FireOptions {
    /// Bam HiFi file with m6A and MSP calls
    #[clap(default_value = "-")]
    pub bam: String,
    /// Output file (bam by default, table if --feats_to_text is used, and bed9 + if --extract is used)
    #[clap(default_value = "-")]
    pub out: String,
    /// Output just FIRE elements in bed9 format
    #[clap(short, long)]
    pub extract: bool,
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
    #[clap(short, long, default_value = "40", env)]
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
    #[clap(short, long, default_value = "85", env)]
    pub min_msp_length_for_positive_fire_call: i64,
    /// Optional path to a model json file.
    /// If not provided ft will use the default model (recommended).
    #[clap(long, env = "FIRE_MODEL")]
    pub model: Option<String>,
    /// Optional path to a FDR table
    #[clap(long, env)]
    pub fdr_table: Option<String>,
    /// Output FIREs features for training in a table format
    #[clap(short, long)]
    pub feats_to_text: bool,
}

#[derive(Args, Debug)]
pub struct AddNucleosomeOptions {
    /// Bam HiFi file with m6A calls
    #[clap(default_value = "-")]
    pub bam: String,
    /// Output bam file with nucleosome calls
    #[clap(default_value = "-")]
    pub out: String,
    /// Minium nucleosome length
    #[clap(short, long, default_value = NUC_LEN)]
    pub nucleosome_length: i64,
    /// Minium nucleosome length when combining over a single m6A
    #[clap(short, long, default_value = COMBO_NUC_LEN)]
    pub combined_nucleosome_length: i64,
    /// Minium distance needed to add to an already existing nuc by crossing an m6a
    #[clap(short, long, default_value = MIN_DIST_ADDED)]
    pub min_distance_added: i64,
    /// Minimum distance from the end of a fiber to call a nucleosome or MSP
    #[clap(short, long, default_value = DIST_FROM_END)]
    pub distance_from_end: i64,
    /// Most m6A events we can skip over to get to the nucleosome length when using D-segment algorithm. 2 is often a good value, negative values disable D-segment for the simple caller.
    #[clap(short, long, default_value = ALLOWED_SKIPS, hide = true)]
    pub allowed_m6a_skips: i64,
    /// Minium score in the ML tag to use in predicting nucleosomes
    #[clap(long, default_value = MIN_ML_SCORE)]
    pub min_ml_score: u8,
}
impl AddNucleosomeOptions {
    pub fn default_value() -> Self {
        Self {
            bam: "-".to_string(),
            out: "-".to_string(),
            nucleosome_length: NUC_LEN.parse().unwrap(),
            combined_nucleosome_length: COMBO_NUC_LEN.parse().unwrap(),
            min_distance_added: MIN_DIST_ADDED.parse().unwrap(),
            distance_from_end: DIST_FROM_END.parse().unwrap(),
            allowed_m6a_skips: ALLOWED_SKIPS.parse().unwrap(),
            min_ml_score: MIN_ML_SCORE.parse().unwrap(),
        }
    }
}

#[derive(Args)]
pub struct DecoratorOptions {
    /// Bam HiFi file with m6A calls
    #[clap(default_value = "-")]
    pub bam: String,
    /// Output path for bed12 file to be decorated
    #[clap(short, long)]
    pub bed12: String,
    /// Output path for decorator bed file
    #[clap(short, long, default_value = "-")]
    pub decorator: String,
    /// Minium score in the ML tag to include in the output
    #[clap(short, long, default_value = MIN_ML_SCORE)]
    pub min_ml_score: u8,
}

#[derive(Args, Debug, PartialEq, Eq)]
pub struct FootprintOptions {
    /// Indexed and aligned bam file with m6A and MSP calls
    pub bam: String,
    /// BED file with the regions to footprint. Should all contain the same motif with proper strand information, and ideally be ChIP-seq peaks.
    pub bed: String,
    /// yaml describing the modules of the footprint
    pub yaml: String,
    /// Output bam
    #[clap(short, long, default_value = "-")]
    pub out: String,
}

#[derive(Args, Debug, PartialEq, Eq)]
pub struct CenterOptions {
    /// Aligned Fiber-seq bam file.
    pub bam: String,
    /// Bed file on which to center fiberseq reads. Data is adjusted to the start position of the bed file and corrected for strand if the strand is indicated in the 6th column of the bed file. The 4th column will also be checked for the strand but only after the 6th is.        
    /// If you include strand information in the 4th (or 6th) column it will orient data accordingly and use the end position of bed record instead of the start if on the minus strand. This means that profiles of motifs in both the forward and minus orientation will align to the same central position.
    pub bed: String,
    /// Minium score in the ML tag to include in the output
    #[clap(short, long, default_value = MIN_ML_SCORE)]
    pub min_ml_score: u8,
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

#[derive(Args, Debug, PartialEq, Eq)]
pub struct ExtractOptions {
    /// Input fiberseq bam file. If no path is provided extract will read bam data from stdin.
    #[arg(default_value = "-", value_hint = ValueHint::AnyPath)]
    pub bam: String,
    /// Report in reference sequence coordinates
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
    /// Minium score in the ML tag to include in the output
    #[clap(short, long, default_value = MIN_ML_SCORE)]
    pub min_ml_score: u8,
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
    /// Simplify output by remove fiber sequence
    #[clap(short, long, help_heading = "All-Format-Options")]
    pub simplify: bool,
}

#[derive(Args, Debug, PartialEq, Eq)]
pub struct ClearKineticsOptions {
    /// Bam HiFi file with kinetics
    #[clap(default_value = "-")]
    pub bam: String,
    /// Output bam file without hifi kinetics
    #[clap(default_value = "-")]
    pub out: String,
}

#[derive(Args, Debug, PartialEq, Eq)]
pub struct StripBasemodsOptions {
    /// Bam HiFi file with base mods
    #[clap(default_value = "-")]
    pub bam: String,
    /// Output bam file
    #[clap(default_value = "-")]
    pub out: String,
    #[clap(short, long, default_value = "m6A",  value_parser(["m6A","6mA", "5mC","CpG"]))]
    /// base modification to strip out of the bam file
    pub basemod: String,
}

#[derive(Args, Debug, PartialEq, Eq)]
pub struct CompletionOptions {
    /// If provided, outputs the completion file for given shell
    #[arg(value_enum)]
    pub shell: Shell,
}

#[derive(Args, Debug, PartialEq, Eq)]
pub struct PileupOptions {
    /// Indexed bam HiFi file with m6A and MSP calls
    pub bam: String,
    /// Region string to make a pileup of.
    /// If not provided will make a pileup of the whole genome
    #[clap(default_value = None)]
    pub rgn: Option<String>,
    /// Output file
    #[clap(short, long, default_value = "-")]
    pub out: String,
    /// Keep zero coverage regions
    #[clap(short, long)]
    pub keep_zeros: bool,
}

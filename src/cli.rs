use anstyle;
use clap::{Args, Command, CommandFactory, Parser, Subcommand};
use clap_complete::{generate, Generator, Shell};
use std::{fmt::Debug, io};

// Reference the modules for the subcommands
mod center_opts;
mod clear_kinetics_opts;
mod ddda_to_m6a_opts;
mod decorator_opts;
mod extract_opts;
mod fiber_hmm;
mod fire_opts;
mod footprint_opts;
mod nucleosome_opts;
mod pg_inject_opts;
mod pg_lift_opts;
mod pg_pansn_opts;
mod pileup_opts;
mod predict_opts;
mod qc_opts;
mod strip_basemods_opts;
mod validate_opts;

// include the subcommand modules as top level functions and structs in the cli module
pub use center_opts::*;
pub use clear_kinetics_opts::*;
pub use ddda_to_m6a_opts::*;
pub use decorator_opts::*;
pub use extract_opts::*;
pub use fiber_hmm::*;
pub use fire_opts::*;
pub use footprint_opts::*;
pub use nucleosome_opts::*;
pub use pg_inject_opts::*;
pub use pg_lift_opts::*;
pub use pg_pansn_opts::*;
pub use pileup_opts::*;
pub use predict_opts::*;
pub use qc_opts::*;
pub use strip_basemods_opts::*;
pub use validate_opts::ValidateOptions;

//
// The main CLI structure
//
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
#[command(version = &**crate::FULL_VERSION)]
#[command(styles=get_styles())]
pub struct Cli {
    #[clap(flatten)]
    pub global: GlobalOpts,
    /// Subcommands for fibertools-rs
    #[clap(subcommand)]
    pub command: Option<Commands>,
}

//
// Global options available to all subcommands
//
#[derive(Debug, Args, PartialEq, Eq)]
pub struct GlobalOpts {
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
}

impl std::default::Default for GlobalOpts {
    fn default() -> Self {
        Self {
            threads: 8,
            verbose: 0,
            quiet: false,
        }
    }
}

//
// This structure contains all the subcommands and their help descriptions.
//
#[derive(Subcommand, Debug)]
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
    /// See https://fiberseq.github.io/fibertools/extracting/extract.html for a description of the outputs.
    #[clap(visible_aliases = &["ex", "e"])]
    Extract(ExtractOptions),
    /// This command centers fiberseq data around given reference positions. This is useful for making aggregate m6A and CpG observations, as well as visualization of SVs.
    ///
    ///  See https://fiberseq.github.io/fibertools/extracting/center.html for a description of the output.
    #[clap(visible_aliases = &["c", "ct"])]
    Center(CenterOptions),
    /// Infer footprints from fiberseq data
    Footprint(FootprintOptions),
    /// Collect QC metrics from a fiberseq bam file
    Qc(QcOpts),
    /// Make decorated bed files for fiberseq data
    TrackDecorators(DecoratorOptions),
    /// Make a pileup track of Fiber-seq features from a FIRE bam
    Pileup(PileupOptions),
    /// Remove HiFi kinetics tags from the input bam file
    ClearKinetics(ClearKineticsOptions),
    /// Strip out select base modifications
    StripBasemods(StripBasemodsOptions),
    /// Convert a DddA BAM file to pseudo m6A BAM file
    DddaToM6a(DddaToM6aOptions),
    /// Apply FiberHMM to a bam file
    FiberHmm(FiberHmmOptions),
    /// Validate a Fiber-seq BAM file for m6A, nucleosome, and optionally FIRE calls
    Validate(ValidateOptions),
    /// Create a mock BAM file from a reference FASTA with perfectly aligned sequences
    #[clap(name = "pg-inject")]
    PgInject(PgInjectOptions),
    /// Lift annotations through a pangenome graph from source to target coordinates
    #[clap(name = "pg-lift")]
    PgLift(PgLiftOptions),
    /// Add or strip panSN-spec prefixes from BAM contig names
    #[clap(name = "pg-pansn")]
    PgPansn(PgPansnOptions),
    /// Make command line completions
    #[clap(hide = true)]
    Completions(CompletionOptions),
    /// Make a man page for fibertools-rs
    ///
    /// Writes file for `man` to stdout.
    #[clap(hide = true)]
    Man {},
}

//
// CLI utility functions
//

/// This function is used to generate the styles for the CLI help messages.
fn get_styles() -> clap::builder::Styles {
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
pub struct CompletionOptions {
    /// If provided, outputs the completion file for given shell
    #[arg(value_enum)]
    pub shell: Shell,
}

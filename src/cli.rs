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
    #[clap(short, long, default_value_t = 8)]
    pub threads: usize,

    /// Logging level [-v: Info, -vv: Debug, -vvv: Trace].
    #[clap(short, long, parse(from_occurrences), help_heading = "DEBUG")]
    pub verbose: usize,

    #[clap(subcommand)]
    pub command: Option<Commands>,
}

///
/// This structure contains all the subcommands for rustybam and their help descriptions.
///
/// Because of naming conventions for rust enums the commands names have
/// different capitalization than on the command line.
/// For example, the `Liftover` enum is invoked using `rustybam liftover`
/// and the `TrimPaf` command with `rustybam trim-paf`.
///
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Trim PAF records that overlap in query sequence to find and optimal splitting point using dynamic programing.
    ///
    /// Note, this can be combined with `rb invert` to also trim the target sequence.
    ///
    /// This idea is to mimic some of the trimming that happens in PAV to improve breakpoint detection. Starts with the largest overlap and iterates until no query overlaps remain.
    #[clap(visible_aliases = &["extract", "e"])]
    Extract {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        bam: String,
        /// Remove contained alignments as well as overlaps.
        #[clap(short, long)]
        reference: bool,
    },
}

pub fn make_cli_parse() -> Cli {
    Cli::parse()
}

pub fn make_cli_app() -> clap::Command<'static> {
    Cli::command()
}

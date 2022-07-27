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
        /// Output path for m6a bed12.
        #[clap(short, long)]
        m6a: Option<String>,
        /// Output path for CpG (primrose) bed12.
        #[clap(short, long)]
        cpg: Option<String>,
        /// Output path for methylation sensitive patch (msp) bed12.
        #[clap(short, long)]
        msp: Option<String>,
        /// Output path for nucleosome bed12.
        #[clap(short, long)]
        nuc: Option<String>,
    },
}

pub fn make_cli_parse() -> Cli {
    Cli::parse()
}

pub fn make_cli_app() -> clap::Command<'static> {
    Cli::command()
}

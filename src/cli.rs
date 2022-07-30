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

    /// Logging level [-v: Debug, -vv: Trace].
    #[clap(short, long, parse(from_occurrences), help_heading = "DEBUG")]
    pub verbose: usize,
    /// Turn of all logging.
    #[clap(short, long)]
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
        /// Output path for m6a bed12.
        #[clap(long)]
        m6a: Option<String>,
        /// Output path for CpG (primrose) bed12.
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
    },
    /// This command centers fiberseq data around given reference positions. This is useful for making aggregate m6a and CpG observations, as well as visualization of SVs.
    Center {
        /// fiberseq bam file, must be aligned and have an index.
        bam: String,
        /// Bed file on which to center fiberseq reads. Data is adjusted to the start position of the bed file. However, if you include strand information in the 4th column it will orient data accordingly and use the end position of bed record instead of the start if on the minus strand. This means that profiles of motifs in both the forward and minus orientation will align to the same central position.
        bed: String,
    },
}

pub fn make_cli_parse() -> Cli {
    Cli::parse()
}

pub fn make_cli_app() -> clap::Command<'static> {
    Cli::command()
}

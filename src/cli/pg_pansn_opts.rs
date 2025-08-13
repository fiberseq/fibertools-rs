use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
pub struct PgPansnOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output BAM file
    #[clap(default_value = "-")]
    pub out: String,
    /// panSN-spec prefix to add to contig names (e.g., "HG002#1#")
    /// If provided, adds the prefix. Mutually exclusive with --strip
    #[clap(short, long, conflicts_with = "strip")]
    pub prefix: Option<String>,
    /// Strip panSN-spec information from contig names
    /// Mutually exclusive with --prefix
    #[clap(short, long, conflicts_with = "prefix")]
    pub strip: bool,
    /// Delimiter character to use when stripping panSN-spec (default: '#')
    #[clap(short, long, default_value = "#", requires = "strip")]
    pub delimiter: char,
    /// Tag to identify haplotype 1 contigs (e.g., "haplotype1")
    /// If both hap1-tag and hap2-tag are provided, reads will be tagged with HP based on contig names
    #[clap(short = '1', long = "hap1-tag")]
    pub hap1_tag: Option<String>,
    /// Tag to identify haplotype 2 contigs (e.g., "haplotype2")
    /// If both hap1-tag and hap2-tag are provided, reads will be tagged with HP based on contig names
    #[clap(short = '2', long = "hap2-tag")]
    pub hap2_tag: Option<String>,
    /// The mapping quality must be greater than or equal to this value for haplotag assignment
    #[clap(short = 'q', long = "min-mapq", default_value = "0")]
    pub min_mapq: u8,
}

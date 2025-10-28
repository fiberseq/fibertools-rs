use crate::utils::input_bam::InputBam;
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug, Clone)]
pub struct PansnParameters {
    /// panSN-spec prefix to add to contig names (e.g., "HG002#1#")
    /// If provided, adds the prefix. Mutually exclusive with --strip
    #[clap(short, long, conflicts_with = "strip")]
    pub prefix: Option<String>,
    /// Strip panSN-spec information from contig names
    /// Mutually exclusive with --prefix
    #[clap(long, conflicts_with = "prefix")]
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
    /// BAM file to copy header from (excluding SQ and HD tags)
    #[clap(short = 'c', long = "copy-header")]
    pub copy_header: Option<String>,
}

impl Default for PansnParameters {
    fn default() -> Self {
        Self {
            prefix: None,
            strip: false,
            delimiter: '#',
            hap1_tag: None,
            hap2_tag: None,
            min_mapq: 0,
            copy_header: None,
        }
    }
}

impl PansnParameters {
    /// Check if any panSN operation is specified in the parameters
    pub fn has_operations(&self) -> bool {
        let has_panspec_operation = self.prefix.is_some() || self.strip;
        let has_haplotag_operation = self.hap1_tag.is_some() && self.hap2_tag.is_some();
        let has_copy_header_operation = self.copy_header.is_some();

        has_panspec_operation || has_haplotag_operation || has_copy_header_operation
    }
}

#[derive(Args, Debug)]
pub struct PgPansnOptions {
    #[clap(flatten)]
    pub input: InputBam,
    /// Output BAM file
    #[clap(default_value = "-")]
    pub out: String,
    /// Additionally write the output BAM header to this file
    #[clap(long)]
    pub header_out: Option<String>,
    #[clap(flatten)]
    pub pansn: PansnParameters,
}

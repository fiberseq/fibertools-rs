use crate::cli;
use crate::fiber::FiberseqRecords;
use crate::utils::bio_io;
use clap::{Args, ValueHint};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::fmt::Debug;

pub static MIN_ML_SCORE: &str = "125";

/// This struct establishes the way a Fiber-seq bam should be filtered when it is read in.
/// This struct is both used as an argument for building FiberSeq records and as a struct that is parsed into cli arguments.
#[derive(Debug, Args, Clone)]
pub struct FiberFilters {
    /// BAM bit flags to filter on, equivalent to `-F` in samtools view
    #[clap(
        global = true,
        short = 'F',
        long = "filter",
        default_value = "0",
        help_heading = "BAM-Options"
    )]
    pub bit_flag: u16,
    /// Filtering expression to use for filtering records
    /// Example: filter to nucleosomes with lengths greater than 150 bp
    ///   -x "len(nuc)>150"
    /// Example: filter to msps with lengths between 30 and 49 bp
    ///   -x "len(msp)=30:50"
    /// Example: combine 2+ filter expressions
    ///   -x "len(nuc)<150,len(msp)=30:50"
    /// Filtering expressions support len() and qual() functions over msp, nuc, m6a, cpg
    #[clap(
        global = true,
        short = 'x',
        long = "ftx",
        alias = "ft-expression",
        help_heading = "BAM-Options"
    )]
    pub filter_expression: Option<String>,
    /// Minium score in the ML tag to use or include in the output
    #[clap(long="ml", alias="min-ml-score", default_value = MIN_ML_SCORE, help_heading = "BAM-Options", env="FT_MIN_ML_SCORE")]
    pub min_ml_score: u8,
    /// Output uncompressed BAM files
    #[clap(help_heading = "BAM-Options", short, long)]
    pub uncompressed: bool,
    /// strip basemods in the first or last X bp of the read
    #[clap(
        global = true,
        long,
        default_value = "0",
        help_heading = "BAM-Options",
        hide = true
    )]
    pub strip_starting_basemods: i64,
}

impl std::default::Default for FiberFilters {
    fn default() -> Self {
        Self {
            bit_flag: 0,
            min_ml_score: MIN_ML_SCORE.parse().unwrap(),
            filter_expression: None,
            uncompressed: false,
            strip_starting_basemods: 0,
        }
    }
}

impl FiberFilters {
    /// This function accepts an iterator over bam records and filters them based on the bit flag.
    pub fn filter_on_bit_flags<'a, I>(
        &'a self,
        records: I,
    ) -> impl Iterator<Item = bam::Record> + 'a
    where
        I: IntoIterator<Item = Result<bam::Record, rust_htslib::errors::Error>> + 'a,
    {
        records
            .into_iter()
            .map(|r| r.expect("htslib is unable to read a record in the input."))
            .filter(|r| {
                // filter by bit flag
                (r.flags() & self.bit_flag) == 0
            })
    }
}

/// This struct is used to parse the input bam file and the filters that should be applied to the bam file.
/// This struct is parsed to create command line arguments and then passed to many functions.
#[derive(Debug, Args)]
pub struct InputBam {
    /// Input BAM file. If no path is provided stdin is used. For m6A prediction, this should be a HiFi bam file with kinetics data. For other commands, this should be a bam file with m6A calls.
    #[clap(default_value = "-", value_hint = ValueHint::AnyPath)]
    pub bam: String,
    #[clap(flatten)]
    pub filters: FiberFilters,
    #[clap(flatten)]
    pub global: cli::GlobalOpts,
    /// by skipping this field it is not parsed as a command line argument
    #[clap(skip)]
    pub header: Option<bam::Header>,
}

impl InputBam {
    pub fn bam_reader(&mut self) -> bam::Reader {
        let mut bam = bio_io::bam_reader(&self.bam);
        bam.set_threads(self.global.threads)
            .expect("unable to set threads for bam reader");
        self.header = Some(bam::Header::from_template(bam.header()));
        bam
    }

    pub fn indexed_bam_reader(&mut self) -> bam::IndexedReader {
        if &self.bam == "-" {
            panic!("Cannot use stdin (\"-\") for indexed bam reading. Please provide a file path for the bam file.");
        }

        let mut bam =
            bam::IndexedReader::from_path(&self.bam).expect("unable to open indexed bam file");
        self.header = Some(bam::Header::from_template(bam.header()));
        bam.set_threads(self.global.threads).unwrap();
        bam
    }

    pub fn fibers<'a>(&self, bam: &'a mut bam::Reader) -> FiberseqRecords<'a> {
        FiberseqRecords::new(bam, self.filters.clone())
    }

    pub fn header_view(&self) -> bam::HeaderView {
        bam::HeaderView::from_header(self.header.as_ref().expect(
            "Input bam must be opened before opening the header or creating a writer with the input bam as a template.",
        ))
    }

    pub fn header(&self) -> bam::Header {
        bam::Header::from_template(&self.header_view())
    }

    pub fn bam_writer(&self, out: &str) -> bam::Writer {
        let header = self.header();
        let program_name = "fibertools-rs";
        let program_id = "ft";
        let program_version = crate::VERSION;
        let mut out = crate::utils::bio_io::program_bam_writer_from_header(
            out,
            header,
            program_name,
            program_id,
            program_version,
        );
        out.set_threads(self.global.threads)
            .expect("unable to set threads for bam writer");
        if self.filters.uncompressed {
            out.set_compression_level(bam::CompressionLevel::Uncompressed)
                .expect("Unable to set compression level to uncompressed");
        }
        out
    }

    /// Add panSN-spec prefix to all contig names in the stored header
    /// This must be called after opening the BAM reader to have an effect
    pub fn add_pansn_prefix(&mut self, pansn_prefix: &str) {
        if let Some(ref mut header) = self.header {
            *header = crate::utils::panspec::add_pan_spec_header(header, pansn_prefix);
        }
    }

    /// Strip panSN-spec information from all contig names in the stored header
    /// using the provided delimiter character (e.g., '#' for panSN format)
    /// This must be called after opening the BAM reader to have an effect
    pub fn strip_pansn_spec(&mut self, delimiter: char) {
        if let Some(ref mut header) = self.header {
            *header = crate::utils::panspec::strip_pan_spec_header(header, &delimiter);
        }
    }
}

impl std::default::Default for InputBam {
    fn default() -> Self {
        Self {
            bam: "-".to_string(),
            filters: FiberFilters::default(),
            global: cli::GlobalOpts::default(),
            header: None,
        }
    }
}

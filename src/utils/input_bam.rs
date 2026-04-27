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
    /// Defaults to 0 (no filtering)
    #[clap(
        global = true,
        short = 'F',
        long = "filter",
        help_heading = "BAM-Options"
    )]
    pub bit_flag: Option<u16>,
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
    /// Convenience: apply the FIRE peak-calling pipeline's fiber-level filters
    /// (`--skip-no-m6a`, `--min-msp 10`, `--min-ave-msp-size 10`). Individual
    /// filter flags still override when both are set. Requires MSP/m6A
    /// annotations on the input BAM, so it is a no-op for commands that run
    /// before those annotations exist.
    #[clap(global = true, long, help_heading = "FIRE-Filter")]
    pub fire_filter: bool,
    /// Drop fibers with no m6A calls. Off by default;
    /// `--fire-filter` turns this on unless explicitly set to `false`.
    /// Use `--skip-no-m6a=false` to override when `--fire-filter` is set.
    #[clap(
        global = true,
        long,
        num_args = 0..=1,
        default_missing_value = "true",
        require_equals = true,
        help_heading = "FIRE-Filter"
    )]
    pub skip_no_m6a: Option<bool>,
    /// Drop fibers with fewer than `N` MSP calls.
    /// Off (0) by default; `--fire-filter` sets this to 10 unless overridden.
    #[clap(global = true, long, env = "MIN_MSP", help_heading = "FIRE-Filter")]
    pub min_msp: Option<usize>,
    /// Drop fibers whose average MSP size is below `N`.
    /// Off (0) by default; `--fire-filter` sets this to 10 unless overridden.
    #[clap(
        global = true,
        long,
        env = "MIN_AVE_MSP_SIZE",
        help_heading = "FIRE-Filter"
    )]
    pub min_ave_msp_size: Option<i64>,
}

impl std::default::Default for FiberFilters {
    fn default() -> Self {
        Self {
            bit_flag: Some(0),
            min_ml_score: MIN_ML_SCORE.parse().unwrap(),
            filter_expression: None,
            uncompressed: false,
            strip_starting_basemods: 0,
            fire_filter: false,
            skip_no_m6a: None,
            min_msp: None,
            min_ave_msp_size: None,
        }
    }
}

impl FiberFilters {
    /// Get the bit flag value, using a default if not explicitly set
    pub fn get_bit_flag(&self) -> u16 {
        self.bit_flag.unwrap_or(0)
    }

    /// Resolved `--skip-no-m6a`, using `--fire-filter` as the fallback.
    pub fn resolved_skip_no_m6a(&self) -> bool {
        self.skip_no_m6a.unwrap_or(self.fire_filter)
    }

    /// Resolved `--min-msp`, using `--fire-filter` (10) as the fallback.
    pub fn resolved_min_msp(&self) -> usize {
        self.min_msp
            .unwrap_or(if self.fire_filter { 10 } else { 0 })
    }

    /// Resolved `--min-ave-msp-size`, using `--fire-filter` (10) as the fallback.
    pub fn resolved_min_ave_msp_size(&self) -> i64 {
        self.min_ave_msp_size
            .unwrap_or(if self.fire_filter { 10 } else { 0 })
    }

    /// True if any FIRE fiber-level filter is active.
    pub fn fire_filter_active(&self) -> bool {
        self.resolved_skip_no_m6a()
            || self.resolved_min_msp() > 0
            || self.resolved_min_ave_msp_size() > 0
    }

    /// True if `rec` passes the FIRE fiber-level filters (skip_no_m6a,
    /// min_msp, min_ave_msp_size). Called by `FiberseqRecords::next` so
    /// every downstream consumer sees only filtered fibers.
    ///
    /// Edge cases worth preserving:
    /// - Fibers with zero MSPs are always rejected once any filter is
    ///   active. The check also guards the divide-by-zero in the average
    ///   MSP size below, so don't drop it when refactoring.
    /// - The no-m6a rejection is gated on `resolved_skip_no_m6a()` so
    ///   `--skip-no-m6a=false` actually disables it (e.g. when combined
    ///   with `--fire-filter`).
    pub fn passes_fire_filter(&self, rec: &crate::fiber::FiberseqData) -> bool {
        if !self.fire_filter_active() {
            return true;
        }
        let n_msps = rec.msp.annotations.len();
        if n_msps == 0 {
            return false;
        }
        if self.resolved_skip_no_m6a() && rec.m6a.annotations.is_empty() {
            return false;
        }
        if n_msps < self.resolved_min_msp() {
            return false;
        }
        let ave_msp_size = rec.msp.lengths().iter().sum::<i64>() / n_msps as i64;
        if ave_msp_size < self.resolved_min_ave_msp_size() {
            return false;
        }
        true
    }

    /// This function accepts an iterator over bam records and filters them based on the bit flag.
    pub fn filter_on_bit_flags<'a, I>(
        &'a self,
        records: I,
    ) -> impl Iterator<Item = bam::Record> + 'a
    where
        I: IntoIterator<Item = Result<bam::Record, rust_htslib::errors::Error>> + 'a,
    {
        let bit_flag = self.get_bit_flag();
        records
            .into_iter()
            .map(|r| r.expect("htslib is unable to read a record in the input."))
            .filter(move |r| {
                // filter by bit flag
                // `move` is needed to capture bit_flag value in the closure
                (r.flags() & bit_flag) == 0
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

    /// Fetch fibers from a specific region with filters applied
    /// Returns an iterator of FiberseqData records from the specified region
    ///
    /// # Arguments
    /// * `bam` - Mutable reference to an IndexedReader
    /// * `chrom` - Chromosome/contig name
    /// * `start` - Optional start position (0-based)
    /// * `end` - Optional end position (0-based, exclusive)
    /// ```
    pub fn fetch_fibers<'a>(
        &'a self,
        bam: &'a mut bam::IndexedReader,
        chrom: &str,
        start: Option<i64>,
        end: Option<i64>,
    ) -> Result<FiberseqRecords<'a, bam::IndexedReader>, rust_htslib::errors::Error> {
        // Fetch the region
        match (start, end) {
            (Some(s), Some(e)) => bam.fetch((chrom, s, e))?,
            (None, None) => bam.fetch(chrom.as_bytes())?,
            _ => panic!("Both start and end must be specified, or neither"),
        }

        // Create FiberseqRecords iterator from the fetched records
        let records = bam.records();
        let header = self.header_view();
        let fiber_iter = FiberseqRecords::from_rec_iterator(records, header, self.filters.clone());

        Ok(fiber_iter)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fiber::FiberseqData;
    use crate::utils::bamannotations::{FiberAnnotation, FiberAnnotations};
    use crate::utils::basemods::BaseMods;

    /// Build a minimal `FiberseqData` shaped only for `passes_fire_filter`.
    /// `msp_lengths` populates `rec.msp.annotations` (each entry contributes
    /// to count and average size). `m6a_count` controls whether `m6a`
    /// annotations are empty or present (only emptiness matters for the
    /// filter).
    fn make_fsd(msp_lengths: &[i64], m6a_count: usize) -> FiberseqData {
        let msp_anns: Vec<FiberAnnotation> = msp_lengths
            .iter()
            .map(|&len| FiberAnnotation {
                start: 0,
                end: len,
                length: len,
                qual: 0,
                reference_start: None,
                reference_end: None,
                reference_length: None,
                extra_columns: None,
            })
            .collect();
        let m6a_anns: Vec<FiberAnnotation> = (0..m6a_count)
            .map(|i| FiberAnnotation {
                start: i as i64,
                end: i as i64 + 1,
                length: 1,
                qual: 0,
                reference_start: None,
                reference_end: None,
                reference_length: None,
                extra_columns: None,
            })
            .collect();
        FiberseqData {
            record: rust_htslib::bam::Record::new(),
            msp: FiberAnnotations::from_annotations(msp_anns, 1000, false),
            nuc: FiberAnnotations::from_annotations(vec![], 1000, false),
            m6a: FiberAnnotations::from_annotations(m6a_anns, 1000, false),
            cpg: FiberAnnotations::from_annotations(vec![], 1000, false),
            base_mods: BaseMods { base_mods: vec![] },
            ec: 0.0,
            target_name: ".".to_string(),
            rg: ".".to_string(),
            center_position: None,
        }
    }

    fn filters() -> FiberFilters {
        FiberFilters::default()
    }

    #[test]
    fn passes_when_no_filter_active() {
        // Default config: nothing rejected, even fibers with no m6a / no msps.
        let f = filters();
        assert!(!f.fire_filter_active());
        assert!(f.passes_fire_filter(&make_fsd(&[], 0)));
        assert!(f.passes_fire_filter(&make_fsd(&[100, 100, 100], 5)));
    }

    #[test]
    fn skip_no_m6a_rejects_only_no_m6a_fibers() {
        let mut f = filters();
        f.skip_no_m6a = Some(true);
        assert!(f.fire_filter_active());
        // No m6a → reject (and no-msp also rejected as a side effect).
        assert!(!f.passes_fire_filter(&make_fsd(&[100], 0)));
        // Has m6a, even with zero msps, still rejected by the empty-msp guard.
        assert!(!f.passes_fire_filter(&make_fsd(&[], 3)));
        // Has both → passes (other thresholds default to 0).
        assert!(f.passes_fire_filter(&make_fsd(&[5], 1)));
    }

    #[test]
    fn min_msp_rejects_fibers_with_too_few_msps() {
        let mut f = filters();
        f.min_msp = Some(3);
        assert!(f.fire_filter_active());
        assert!(!f.passes_fire_filter(&make_fsd(&[100, 100], 5)));
        assert!(f.passes_fire_filter(&make_fsd(&[100, 100, 100], 5)));
    }

    #[test]
    fn min_ave_msp_size_rejects_low_average() {
        let mut f = filters();
        f.min_ave_msp_size = Some(50);
        // Average = 30 → reject.
        assert!(!f.passes_fire_filter(&make_fsd(&[10, 20, 60], 5)));
        // Average = 60 → pass.
        assert!(f.passes_fire_filter(&make_fsd(&[40, 60, 80], 5)));
    }

    #[test]
    fn fire_filter_combo_applies_all_three_defaults() {
        // `--fire-filter` alone should imply skip_no_m6a + min_msp=10 + min_ave_msp_size=10.
        let mut f = filters();
        f.fire_filter = true;
        assert!(f.resolved_skip_no_m6a());
        assert_eq!(f.resolved_min_msp(), 10);
        assert_eq!(f.resolved_min_ave_msp_size(), 10);
        // <10 msps → reject.
        let nine_long_msps: Vec<i64> = vec![100; 9];
        assert!(!f.passes_fire_filter(&make_fsd(&nine_long_msps, 5)));
        // 10 msps but ave_size = 5 < 10 → reject.
        let ten_short_msps: Vec<i64> = vec![5; 10];
        assert!(!f.passes_fire_filter(&make_fsd(&ten_short_msps, 5)));
        // 10 msps, ave_size = 100, has m6a → pass.
        let ten_long_msps: Vec<i64> = vec![100; 10];
        assert!(f.passes_fire_filter(&make_fsd(&ten_long_msps, 5)));
        // Same fiber but no m6a → reject (skip_no_m6a is implied).
        assert!(!f.passes_fire_filter(&make_fsd(&ten_long_msps, 0)));
    }

    #[test]
    fn explicit_flag_overrides_fire_filter_default() {
        // `--fire-filter --min-msp=5` should use 5, not 10.
        let mut f = filters();
        f.fire_filter = true;
        f.min_msp = Some(5);
        assert_eq!(f.resolved_min_msp(), 5);
        // Other two still take fire-filter defaults.
        assert!(f.resolved_skip_no_m6a());
        assert_eq!(f.resolved_min_ave_msp_size(), 10);
        // 5 msps would have failed under the 10 default, but passes under 5.
        let five_long_msps: Vec<i64> = vec![100; 5];
        assert!(f.passes_fire_filter(&make_fsd(&five_long_msps, 5)));
    }

    #[test]
    fn explicit_skip_no_m6a_false_overrides_fire_filter_default() {
        // `--fire-filter --skip-no-m6a=false` keeps the size/count thresholds
        // but should turn the no-m6a guard off. The empty-msp guard still
        // applies (it's a divide-by-zero protection), so we use a record
        // with msps but no m6a.
        let mut f = filters();
        f.fire_filter = true;
        f.skip_no_m6a = Some(false);
        assert!(!f.resolved_skip_no_m6a());
        let ten_long_msps: Vec<i64> = vec![100; 10];
        assert!(f.passes_fire_filter(&make_fsd(&ten_long_msps, 0)));
    }
}

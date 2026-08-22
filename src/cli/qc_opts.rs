use crate::utils::input_bam::{CallableFibers, InputBam};
use clap::Args;
use std::fmt::Debug;

#[derive(Args, Debug)]
#[command(after_long_help = r#"Output schema:
  A tab-separated table with the columns "statistic	value	count	count_filtered".

  The "count" column counts every read in the stream. The stream is already
  smaller if you use -F, --ml, --strip-starting-basemods, or -x. The
  "count_filtered" column counts only Callable reads. The fiberseq_callable
  rows show the three states: Callable, NotCallable, and Untagged.
  Untagged reads never enter the filtered column.

  Six statistics measure bases: fiber_length, phased_bp, m6a_count,
  cpg_count, m6a_ratio, and read_length_per_nuc. For these, the filtered
  column uses only the callable span of each read. Their filtered keys come
  from the span, so a row can have a zero in either column, and
  "count_filtered <= count" does not hold for these rows. The filtered
  m6a_ratio is a little higher than the unfiltered one, because the span
  starts and ends at m6A calls.

  The callability minimums (--min-msp, --min-ave-msp-size, defaults 10
  and 10) set which reads are Callable and so select the filtered
  column. They never drop reads on their own. --callable-fibers
  (alias --fire-filter) restricts the stream to callable fibers, which
  shrinks the count column too. With
  --acf, each m6a_acf row carries the filtered ACF in the 4th field (NA
  when no sampled read passes), and the acf_reads row shows the sample size
  and how many of those reads pass."#)]
pub struct QcOpts {
    #[clap(flatten)]
    pub input: InputBam<CallableFibers>,
    /// Output text file with QC metrics. See the schema notes at the
    /// bottom of --help.
    #[clap(default_value = "-")]
    pub out: String,
    /// Calculate the auto-correlation function of the m6A marks in the fiber-seq data.
    #[clap(long)]
    pub acf: bool,
    /// maximum lag for the ACF calculation
    #[clap(long, default_value = "250")]
    pub acf_max_lag: usize,
    /// Minimum number of m6A marks to use a read in the ACF calculation
    #[clap(long, default_value = "100")]
    pub acf_min_m6a: usize,
    /// maximum number of reads to use in the ACF calculation
    #[clap(long, default_value = "10000")]
    pub acf_max_reads: usize,
    /// After sampling the first "acf-max-reads" randomly sample one of every "acf-sample-rate" reads and replace one of the previous reads at random.
    #[clap(long, default_value = "100")]
    pub acf_sample_rate: f32,
    /// In the output include a measure of the number of m6A events per MSPs of a given size.
    /// The output format is: "m6a_per_msp_size\t{m6A count},{MSP size},{is a FIRE}\t{count}\t{count_filtered}"
    /// e.g. "m6a_per_msp_size\t35,100,false\t100\t98"
    #[clap(short, long)]
    pub m6a_per_msp: bool,
    /// Only process the first "n" reads in the input bam file.
    #[clap(long)]
    pub n_reads: Option<usize>,
    /// Seed for the ACF read sampler, so --acf output is reproducible by
    /// default. Pass --random-seed to sample from entropy instead.
    #[clap(long, default_value = "42", conflicts_with = "random_seed")]
    pub seed: u64,
    /// Use a random seed for the ACF read sampler instead of --seed.
    #[clap(long)]
    pub random_seed: bool,
    /// A read enters the count_filtered column only if its callable SPAN
    /// is at least this many bases. The span is the fiberseq_callable
    /// range, not the read length: a long read with a short callable
    /// span fails this filter. Applies on top of the callable state.
    #[clap(long, value_parser = clap::value_parser!(i64).range(0..))]
    pub filtered_min_callable_length: Option<i64>,
}

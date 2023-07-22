/// Add and remove base modifications from a bam record
pub mod basemods;
/// Center fiberseq information around a reference position
pub mod center;
#[cfg(feature = "cli")]
/// Command line interface for fibertools-rs.
pub mod cli;
/// Extract fiberseq data into plain text formats
pub mod extract;
/// Add nucleosomes to a bam file
pub mod nucleosomes;
#[cfg(feature = "predict")]
/// m6A prediction
pub mod predict_m6a;
/// Remove base modifications from a bam record
pub mod strip_basemods;

use anyhow::Result;
use bio_io::*;
use indicatif::{style, ProgressBar};
#[cfg(feature = "predict")]
use itertools::Itertools;
use rust_htslib::{bam, bam::Read};
use std::env;
use std::io::Write;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const GIT_HASH: &str = env!("CARGO_GIT_HASH");
pub const LONG_VERSION: &str = env!("CARGO_LONG_VERSION");
const PROGRESS_STYLE: &str =
    "[{elapsed_precise:.yellow}] {bar:50.cyan/blue} {human_pos:>5.cyan}/{human_len:.blue} {percent:>3.green}% {per_sec:<10.cyan}";

/// unzip a vector of tuples
pub fn unzip_to_vectors<T, U>(vec: Vec<(T, U)>) -> (Vec<T>, Vec<U>) {
    vec.into_iter().unzip()
}

/// join a vector with commas
pub fn join_by_str<'a, I, Z>(vals: I, sep: &str) -> String
where
    I: IntoIterator<Item = Z>,
    Z: ToString + 'a,
{
    vals.into_iter().map(|v| v.to_string() + sep).collect()
}

pub struct FiberOut {
    pub m6a: Option<Box<dyn Write>>,
    pub cpg: Option<Box<dyn Write>>,
    pub msp: Option<Box<dyn Write>>,
    pub nuc: Option<Box<dyn Write>>,
    pub all: Option<Box<dyn Write>>,
    pub reference: bool,
    pub simplify: bool,
    pub quality: bool,
    pub min_ml_score: u8,
    pub full_float: bool,
}

impl FiberOut {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        m6a: &Option<String>,
        cpg: &Option<String>,
        msp: &Option<String>,
        nuc: &Option<String>,
        all: &Option<String>,
        reference: bool,
        simplify: bool,
        quality: bool,
        min_ml_score: u8,
        full_float: bool,
    ) -> Result<Self> {
        let m6a = match m6a {
            Some(m6a) => Some(writer(m6a)?),
            None => None,
        };
        let cpg = match cpg {
            Some(cpg) => Some(writer(cpg)?),
            None => None,
        };
        let msp = match msp {
            Some(msp) => Some(writer(msp)?),
            None => None,
        };
        let nuc = match nuc {
            Some(nuc) => Some(writer(nuc)?),
            None => None,
        };
        let all = match all {
            Some(all) => Some(writer(all)?),
            None => None,
        };
        // set to zero
        let mut min_ml_score = min_ml_score;
        if full_float {
            min_ml_score = 0;
        }

        Ok(FiberOut {
            m6a,
            cpg,
            msp,
            nuc,
            all,
            reference,
            simplify,
            quality,
            min_ml_score,
            full_float,
        })
    }
}

/// clear kinetics from a hifi bam
pub fn clear_kinetics(bam: &mut bam::Reader, out: &mut bam::Writer) {
    let bar = ProgressBar::new(1);
    let style_str="[Clearing Kinetics] [Elapsed {elapsed:.yellow}] [Reads processed {human_pos:>5.cyan}] (reads/s {per_sec:>5.green})";
    let style = style::ProgressStyle::with_template(style_str)
        .unwrap()
        .progress_chars("##-");
    bar.set_style(style);
    for rec in bam.records() {
        let mut record = rec.unwrap();
        record.remove_aux(b"fp").unwrap_or(());
        record.remove_aux(b"fi").unwrap_or(());
        record.remove_aux(b"rp").unwrap_or(());
        record.remove_aux(b"ri").unwrap_or(());
        out.write(&record).unwrap();
        bar.inc(1);
    }
    bar.finish();
}

/// Write to a bam file.
pub fn bam_writer(out: &str, template_bam: &bam::Reader, threads: usize) -> bam::Writer {
    let program_name = "fibertools-rs";
    let program_id = "ft";
    let program_version = VERSION;
    program_bam_writer(
        out,
        template_bam,
        threads,
        program_name,
        program_id,
        program_version,
    )
}

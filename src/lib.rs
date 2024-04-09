/// Data structure for fiberseq data
pub mod fiber;

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
/// m6A prediction
pub mod predict_m6a;
/// Remove base modifications from a bam record
pub mod strip_basemods;

/// add fire data
pub mod fire;

/// make a fire track from a bam file
pub mod firetrack;

/// add decorators
pub mod decorator;

pub mod footprint;

pub mod bamlift;

pub mod bio_io;

pub mod bamranges;

pub mod m6a_burn;

use anyhow::Result;
use bio_io::*;
use itertools::Itertools;
use rust_htslib::{bam, bam::Read};
use std::env;
use std::io::Write;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const GIT_HASH: &str = env!("CARGO_GIT_HASH");
pub const LONG_VERSION: &str = env!("CARGO_LONG_VERSION");
// if this string (bar)gets too long it displays weird when writing to stdout
const PROGRESS_STYLE: &str =
    "[{elapsed_precise:.yellow}] {bar:>35.cyan/blue} {human_pos:>5.cyan}/{human_len:.blue} {percent:>3.green}% {per_sec:<10.cyan}";

/// COLORS
pub const NUC_COLOR: &str = "169,169,169";
pub const M6A_COLOR: &str = "128,0,128";
pub const CPG_COLOR: &str = "139,69,19";
pub const LINKER_COLOR: &str = "147,112,219";
pub const FIRE_COLORS: [(f32, &str); 9] = [
    (1.0, "139,0,0"),
    (2.0, "175,0,0"),
    (3.0, "200,0,0"),
    (4.0, "225,0,0"),
    (5.0, "255,0,0"),
    (10.0, "255,140,0"),
    (25.0, "225,225,0"),
    (100.0, LINKER_COLOR),
    (200.0, NUC_COLOR),
];

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

/// join a vector with commas
pub fn join_by_str_option(vals: &[Option<i64>], sep: &str) -> String {
    vals.iter()
        .map(|v| match v {
            Some(v) => v.to_string(),
            None => String::from("NA"),
        })
        .map(|v| v + sep)
        .collect()
}

/// join a vector with commas
pub fn join_by_str_option_can_skip(vals: &[Option<i64>], sep: &str, skip_none: bool) -> String {
    vals.iter()
        .map(|v| match v {
            Some(v) => v.to_string(),
            None => {
                if skip_none {
                    String::from("")
                } else {
                    String::from("NA")
                }
            }
        })
        .filter(|v| !v.is_empty())
        .map(|v| v + sep)
        .collect()
}

/// clear kinetics from a hifi bam
pub fn clear_kinetics(bam: &mut bam::Reader, out: &mut bam::Writer) {
    let bar = bio_io::no_length_progress_bar();
    for rec in bam.records() {
        let mut record = rec.unwrap();
        record.remove_aux(b"fp").unwrap_or(());
        record.remove_aux(b"fi").unwrap_or(());
        record.remove_aux(b"rp").unwrap_or(());
        record.remove_aux(b"ri").unwrap_or(());
        out.write(&record).unwrap();
        bar.inc_length(1);
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

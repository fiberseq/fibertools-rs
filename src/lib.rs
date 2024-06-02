/// Add and remove base modifications from a bam record
pub mod basemods;
/// Center fiberseq information around a reference position
pub mod center;
#[cfg(feature = "cli")]
/// Command line interface for fibertools-rs.
pub mod cli;
/// Extract fiberseq data into plain text formats
pub mod extract;
/// Data structure for fiberseq data
pub mod fiber;
/// add fire data
pub mod fire;
/// Add nucleosomes to a bam file
pub mod nucleosomes;
/// make a fire track from a bam file
pub mod pileup;
/// m6A prediction
pub mod predict_m6a;
/// Remove base modifications from a bam record
pub mod strip_basemods;

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
use rust_htslib::bam::FetchDefinition;
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
pub fn clear_kinetics(opts: &mut cli::ClearKineticsOptions) {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    //let mut out = bam_writer(&opts.out, &bam, opts.input.global.threads);

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

pub fn region_parser(rgn: &str) -> (FetchDefinition<'_>, String) {
    if rgn.contains(':') {
        let (chrom, rest) = rgn.split(':').collect_tuple().unwrap();
        let (start, end) = rest.split('-').collect_tuple().unwrap();
        (
            FetchDefinition::RegionString(
                chrom.as_bytes(),
                start.parse().unwrap(),
                end.parse().unwrap(),
            ),
            chrom.to_string(),
        )
    } else {
        (FetchDefinition::String(rgn.as_bytes()), rgn.to_string())
    }
}

/// Data structure for fiberseq data
pub mod fiber;
pub mod m6a_burn;
/// subcommands of fibertools-rs
pub mod subcommands;
pub mod utils;

#[cfg(feature = "cli")]
/// Command line interface for fibertools-rs.
pub mod cli;

use crate::utils::bio_io::*;
use anyhow::Result;
use itertools::Itertools;
use lazy_static::lazy_static;
use rust_htslib::bam::FetchDefinition;
use rust_htslib::{bam, bam::Read};
use std::env;
use std::io::Write;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
lazy_static! {
    pub static ref FULL_VERSION: String = format!(
        "v{}\tgit-details {}",
        env!("CARGO_PKG_VERSION"),
        env!("VERGEN_GIT_DESCRIBE")
    );
}
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

pub fn region_parser(rgn: &str) -> (FetchDefinition<'_>, String) {
    if rgn.contains(':') {
        let (chrom, rest) = rgn.split(':').collect_tuple().unwrap();
        let (start, end) = rest.split('-').collect_tuple().unwrap();
        let st: i64 = start
            .replace(',', "")
            .parse()
            .unwrap_or_else(|_| panic!("Could not parse start of region: {start}"));
        (
            FetchDefinition::RegionString(
                chrom.as_bytes(),
                st,
                end.replace(',', "")
                    .parse()
                    .unwrap_or_else(|_| panic!("Could not parse end of region: {end}")),
            ),
            chrom.to_string(),
        )
    } else {
        (FetchDefinition::String(rgn.as_bytes()), rgn.to_string())
    }
}

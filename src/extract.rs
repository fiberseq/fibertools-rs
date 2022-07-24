use bio::alphabets::dna::revcomp;
use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{
    bam,
    bam::record::Aux,
    bam::{ext::BamRecordExtensions, Read},
};
use std::convert::TryFrom;

pub struct BaseMods {
    pub base: u8,
    pub strand: char,
    pub modification_type: char,
    pub modified_positions: Vec<i64>,
}

impl BaseMods {}

pub fn get_mm_tag(record: &bam::Record) -> Vec<BaseMods> {
    // regex for matching the MM tag
    lazy_static! {
        static ref MM_RE: Regex =
            Regex::new(r"((([ACGTUN])([-+])([a-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
    }
    // Array to store all the different modifications within the MM tag
    let mut rtn = vec![];

    // if there is an MM tag iterate over all the regex matches
    if let Ok(Aux::String(mm_text)) = record.aux(b"MM") {
        for cap in MM_RE.captures_iter(mm_text) {
            let mod_base = cap.get(3).map(|m| m.as_str().as_bytes()[0]).unwrap();
            let mod_strand = cap.get(4).map_or("", |m| m.as_str());
            let modification_type = cap.get(5).map_or("", |m| m.as_str());
            let mod_dists_str = cap.get(6).map_or("", |m| m.as_str());
            // parse the string containing distances between modifications into a vector of i64
            let mod_dists: Vec<i64> = mod_dists_str
                .trim_end_matches(';')
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().unwrap())
                .collect();

            // get forward sequence bases from the bam record
            let forward_bases = if record.is_reverse() {
                revcomp(record.seq().as_bytes())
            } else {
                record.seq().as_bytes()
            };

            // find real positions in the forward sequence
            let mut cur_mod_idx = 0;
            let mut cur_seq_idx = 0;
            let mut dist_from_last_mod_base = 0;
            let mut modified_positions: Vec<i64> = vec![0; mod_dists.len()];
            while cur_seq_idx < forward_bases.len() && cur_mod_idx < mod_dists.len() {
                let cur_base = forward_bases[cur_seq_idx];
                if cur_base == mod_base && dist_from_last_mod_base == mod_dists[cur_mod_idx] {
                    modified_positions[cur_mod_idx] = i64::try_from(cur_seq_idx).unwrap();
                    dist_from_last_mod_base = 0;
                    cur_mod_idx += 1;
                } else if cur_base == mod_base {
                    dist_from_last_mod_base += 1
                }
                cur_seq_idx += 1;
            }
            // assert that we extract the same number of modifications as we have distances
            assert_eq!(cur_mod_idx, mod_dists.len());

            // add to a struct
            let mods = BaseMods {
                base: mod_base,
                strand: mod_strand.chars().next().unwrap(),
                modification_type: modification_type.chars().next().unwrap(),
                modified_positions: modified_positions,
            };
            rtn.push(mods);
        }
    } else {
        log::debug!("No MM tag found");
    }
    return rtn;
}

pub fn extract_from_record(record: &bam::Record, reference: bool) {
    let mods = get_mm_tag(record);
    for moda in mods.iter() {
        if moda.base == b'C' && moda.modification_type == 'm' {
            println!("{:?}", moda.modified_positions);
        }
        // we want to get the bases on the reference sequence when possible
        if reference && !record.is_unmapped() {
            let _aligned_bases = record.aligned_pairs();
        }
    }
}

pub fn extract_contained(bam: &mut bam::Reader, reference: bool) {
    for r in bam.records() {
        let record = r.unwrap();
        extract_from_record(&record, reference);
    }
}

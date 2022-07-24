use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{
    bam,
    bam::record::Aux,
    bam::{ext::BamRecordExtensions, Read},
};
use std::str;

pub struct BaseMods {
    pub base: char,
    pub strand: char,
    pub mod_type: char,
    pub mod_pos: Vec<i64>,
}

pub fn get_mm_tag(record: &bam::Record) -> Vec<BaseMods> {
    lazy_static! {
        static ref MM_RE: Regex =
            Regex::new(r"((([ACGTUN])([-+])([a-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
    }

    let mut rtn = vec![];

    if let Ok(Aux::String(mm_text)) = record.aux(b"MM") {
        //let read_array = array.iter().collect::<Vec<_>>();
        //log::debug!("MM tag: {:?}", mm_text);
        for cap in MM_RE.captures_iter(mm_text) {
            let mod_base = cap
                .get(3)
                .map(|m| m.as_str().chars().next().unwrap())
                .unwrap();
            let mod_strand = cap.get(4).map_or("", |m| m.as_str());
            let mod_type = cap.get(5).map_or("", |m| m.as_str());
            let mod_pos_str = cap.get(6).map_or("", |m| m.as_str());
            let mod_pos: Vec<i64> = mod_pos_str
                .trim_end_matches(';')
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().unwrap())
                .collect();

            log::trace!("Mod: {mod_base}{mod_strand}{mod_type}\nMod pos: {mod_pos:?}\n");
            let mods = BaseMods {
                base: mod_base,
                strand: mod_strand.chars().next().unwrap(),
                mod_type: mod_type.chars().next().unwrap(),
                mod_pos: mod_pos,
            };
            rtn.push(mods);
        }
    } else {
        log::debug!("No MM tag found");
    }
    return rtn;
}

pub fn extract_from_record(record: &bam::Record, reference: bool) {
    if reference && record.is_reverse() {
        println!("{:#?}", str::from_utf8(record.qname()).unwrap());
    }

    get_mm_tag(record);

    if reference && !record.is_unmapped() {
        for [q_pos, r_pos] in record.aligned_pairs() {
            println!("{q_pos}->{r_pos}");
        }
    }
}

pub fn extract_contained(bam: &mut bam::Reader, reference: bool) {
    for r in bam.records() {
        let record = r.unwrap();
        extract_from_record(&record, reference);
    }
}

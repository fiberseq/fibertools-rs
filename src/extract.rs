use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{
    bam,
    bam::record::Aux,
    bam::{ext::BamRecordExtensions, Read},
};
use std::str;

pub fn get_mm_tag(record: &bam::Record) {
    lazy_static! {
        static ref MM_RE: Regex =
            Regex::new(r"(([ACGTUN][-+]([a-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
    }

    if let Ok(Aux::String(mm_text)) = record.aux(b"MM") {
        //let read_array = array.iter().collect::<Vec<_>>();
        //log::debug!("MM tag: {:?}", mm_text);
        for cap in MM_RE.captures_iter(mm_text) {
            let mod_type = cap.get(2).map_or("", |m| m.as_str());
            let mod_pos_str = cap.get(4).map_or("", |m| m.as_str());
            let mod_pos: Vec<i32> = mod_pos_str
                .trim_end_matches(";")
                .split(',')
                .map(|s| s.trim()) // (2)
                .filter(|s| !s.is_empty()) // (3)
                .map(|s| s.parse().unwrap()) // (4)
                .collect();
            log::debug!("Mod type: {mod_type}\nMod pos: {mod_pos:?}\n");
        }
    } else {
        log::debug!("No MM tag found");
    }
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

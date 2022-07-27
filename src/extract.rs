use super::*;
use bio::alphabets::dna::revcomp;
use colored::Colorize;
use lazy_static::lazy_static;
use rayon::{current_num_threads, prelude::*};
use regex::Regex;
use rust_htslib::{bam, bam::record::Aux, bam::Read};
use std::convert::TryFrom;
use std::time::Instant;
pub struct BaseMods {
    pub modified_base: u8,
    pub strand: char,
    pub modification_type: char,
    pub modified_positions: Vec<i64>,
    pub modified_reference_positions: Vec<i64>,
}

impl BaseMods {
    pub fn new(record: &bam::Record) -> Vec<BaseMods> {
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
                let mut mods = BaseMods {
                    modified_base: mod_base,
                    strand: mod_strand.chars().next().unwrap(),
                    modification_type: modification_type.chars().next().unwrap(),
                    modified_positions,
                    modified_reference_positions: vec![],
                };
                // add the reference bases
                mods.add_reference_positions(record);
                rtn.push(mods);
            }
        } else {
            log::debug!("No MM tag found");
        }
        rtn
    }

    pub fn add_reference_positions(&mut self, record: &bam::Record) {
        let positions = positions_on_complimented_sequence(record, &self.modified_positions);
        // get the reference positions
        self.modified_reference_positions = liftover_exact(record, &positions);
    }
}

/// Merge two lists into a sorted list
/// Normal sort is supposed to be very fast on two sorted lists
/// https://doc.rust-lang.org/std/vec/struct.Vec.html#current-implementation-6
pub fn merge_two_lists<T>(left: Vec<T>, right: Vec<T>) -> Vec<T>
where
    T: Ord,
    T: Clone,
{
    let mut x = [left, right].concat();
    x.sort();
    x
}

///```
/// use rust_htslib::{bam, bam::Read};
/// use fibertools_rs::*;
/// use log;
/// use env_logger::{Builder, Target};;
/// Builder::new().target(Target::Stderr).filter(None, log::LevelFilter::Debug).init();
/// let mut bam = bam::Reader::from_path(&".test/aligned.bam").unwrap();
/// for record in bam.records() {
///     let record = record.unwrap();
///     let n_s = extract::get_u32_tag(&record, b"ns");    
///     let n_l = extract::get_u32_tag(&record, b"nl");    
///     let a_s = extract::get_u32_tag(&record, b"as");    
///     let a_l = extract::get_u32_tag(&record, b"al");    
///     log::debug!("{:?}", a_s);
/// }
///```
pub fn get_u32_tag(record: &bam::Record, tag: &[u8; 2]) -> Vec<i64> {
    if let Ok(Aux::ArrayU32(array)) = record.aux(tag) {
        let read_array = array.iter().map(|x| x as i64).collect::<Vec<_>>();
        read_array
    } else {
        vec![]
    }
}

pub struct FiberseqData {
    pub record: bam::Record,
    pub nuc_starts: Vec<i64>,
    pub msp_starts: Vec<i64>,
    pub nuc_length: Vec<i64>,
    pub msp_length: Vec<i64>,
    pub ref_nuc: Vec<(i64, i64)>,
    pub ref_msp: Vec<(i64, i64)>,
    pub base_mods: Vec<BaseMods>,
}

impl FiberseqData {
    pub fn new(record: &bam::Record) -> Self {
        let nuc_starts = get_u32_tag(record, b"ns");
        let msp_starts = get_u32_tag(record, b"as");
        let nuc_length = get_u32_tag(record, b"nl");
        let msp_length = get_u32_tag(record, b"al");
        // range positions
        let ref_nuc = get_closest_reference_range(&nuc_starts, &nuc_length, record);
        let ref_msp = get_closest_reference_range(&msp_starts, &msp_length, record);
        //
        FiberseqData {
            record: record.clone(),
            nuc_starts,
            msp_starts,
            nuc_length,
            msp_length,
            ref_nuc,
            ref_msp,
            base_mods: BaseMods::new(record),
        }
    }

    pub fn from_records(records: &Vec<bam::Record>) -> Vec<Self> {
        records
            .par_iter()
            .map(FiberseqData::new)
            .collect::<Vec<_>>()
    }

    pub fn to_string(&mut self, reference: bool, starts: &[i64], lengths: &[i64]) -> String {
        let ct;
        let start;
        let end;
        let name = std::str::from_utf8(self.record.qname()).unwrap();
        if reference {
            ct = "";
            start = self.record.reference_start();
            end = self.record.reference_end();
        } else {
            ct = name;
            start = 0;
            end = self.record.seq_len() as i64;
        }
        let strand = self.record.is_reverse();
        let score = 0;
        let color = "126,126,126";
        let b_ct = starts.len();
        let b_st: String = starts.iter().map(|&id| id.to_string() + ",").collect();
        let b_ln: String = lengths.iter().map(|&id| id.to_string() + ",").collect();
        format!("{ct}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{start}\t{end}\t{color}\t{b_ct}\t{b_ln}\t{b_st}\n")
    }
}

pub fn process_bam_chunk(
    records: &Vec<bam::Record>,
    so_far: usize,
    _reference: bool,
    m6a: &Option<String>,
    _cpg: &Option<String>,
    _msp: &Option<String>,
    _nuc: &Option<String>,
) {
    let start = Instant::now();
    let _fiber_data = FiberseqData::from_records(records);
    let duration = start.elapsed().as_secs_f64();

    match m6a {
        Some(m6a) => {
            log::info!("Processing {:}", m6a);
        }
        None => {}
    }

    log::info!(
        "Processing {}. {} reads done so far.",
        format!("{:.2?} reads/s", records.len() as f64 / duration)
            .bright_cyan()
            .bold(),
        format!("{:}", so_far + records.len())
            .bright_magenta()
            .bold()
    );
}

pub fn extract_contained(
    bam: &mut bam::Reader,
    reference: bool,
    m6a: &Option<String>,
    cpg: &Option<String>,
    msp: &Option<String>,
    nuc: &Option<String>,
) {
    // process bam in chunks
    // keeps mem pretty low, about 1GB per thread
    let bin_size = current_num_threads() * 2000;
    //
    let mut cur_count = 0;
    let mut cur_vec = vec![];
    let mut proccesed_reads = 0;
    for r in bam.records() {
        let record = r.unwrap();
        cur_vec.push(record);
        cur_count += 1;
        if cur_count == bin_size {
            process_bam_chunk(&cur_vec, proccesed_reads, reference, m6a, cpg, msp, nuc);
            proccesed_reads += cur_vec.len();
            cur_vec.clear();
            cur_count = 0;
        }
    }
    // clear any unprocessed recs not big enough to make a full chunk
    process_bam_chunk(&cur_vec, proccesed_reads, reference, m6a, cpg, msp, nuc);
}

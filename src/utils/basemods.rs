use crate::utils::bamannotations::*;
use crate::utils::bio_io::*;
use bio::alphabets::dna::revcomp;
use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{
    bam,
    bam::record::{Aux, AuxArray},
};
use std::collections::HashMap;

use std::convert::TryFrom;

#[derive(Eq, PartialEq, Debug, PartialOrd, Ord, Clone)]
pub struct BaseMod {
    pub modified_base: u8,
    pub strand: char,
    pub modification_type: char,
    pub ranges: Ranges,
    pub record_is_reverse: bool,
}

impl BaseMod {
    pub fn new(
        record: &bam::Record,
        modified_base: u8,
        strand: char,
        modification_type: char,
        modified_bases_forward: Vec<i64>,
        modified_probabilities_forward: Vec<u8>,
    ) -> Self {
        let tmp = modified_bases_forward.clone();
        let mut ranges = Ranges::new(record, modified_bases_forward, None, None);
        ranges.set_qual(modified_probabilities_forward);
        let record_is_reverse = record.is_reverse();
        assert_eq!(tmp, ranges.forward_starts(), "forward starts not equal");
        Self {
            modified_base,
            strand,
            modification_type,
            ranges,
            record_is_reverse,
        }
    }

    pub fn is_m6a(&self) -> bool {
        self.modification_type == 'a'
    }

    pub fn is_cpg(&self) -> bool {
        self.modification_type == 'm'
    }

    pub fn filter_at_read_ends(&mut self, n_strip: i64) {
        if n_strip <= 0 {
            return;
        }
        self.ranges.filter_starts_at_read_ends(n_strip);
    }
}

#[derive(Eq, PartialEq, Debug, Clone)]
pub struct BaseMods {
    pub base_mods: Vec<BaseMod>,
}

impl BaseMods {
    pub fn new(record: &bam::Record, min_ml_score: u8) -> BaseMods {
        // my basemod parser is ~25% faster than rust_htslib's
        BaseMods::my_mm_ml_parser(record, min_ml_score)
    }

    pub fn my_mm_ml_parser(record: &bam::Record, min_ml_score: u8) -> BaseMods {
        // regex for matching the MM tag
        lazy_static! {
            // MM:Z:([ACGTUN][-+]([A-Za-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*
            static ref MM_RE: Regex =
                Regex::new(r"((([ACGTUN])([-+])([A-Za-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
        }
        // Array to store all the different modifications within the MM tag
        let mut rtn = vec![];

        let ml_tag = get_u8_tag(record, b"ML");

        let mut num_mods_seen = 0;

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
                    convert_seq_uppercase(revcomp(record.seq().as_bytes()))
                } else {
                    convert_seq_uppercase(record.seq().as_bytes())
                };
                log::trace!(
                    "mod_base: {}, mod_strand: {}, modification_type: {}, mod_dists: {:?}",
                    mod_base as char,
                    mod_strand,
                    modification_type,
                    mod_dists
                );
                // find real positions in the forward sequence
                let mut cur_mod_idx = 0;
                let mut cur_seq_idx = 0;
                let mut dist_from_last_mod_base = 0;
                let mut unfiltered_modified_positions: Vec<i64> = vec![0; mod_dists.len()];
                while cur_seq_idx < forward_bases.len() && cur_mod_idx < mod_dists.len() {
                    let cur_base = forward_bases[cur_seq_idx];
                    if (cur_base == mod_base || mod_base == b'N')
                        && dist_from_last_mod_base == mod_dists[cur_mod_idx]
                    {
                        unfiltered_modified_positions[cur_mod_idx] =
                            i64::try_from(cur_seq_idx).unwrap();
                        dist_from_last_mod_base = 0;
                        cur_mod_idx += 1;
                    } else if cur_base == mod_base {
                        dist_from_last_mod_base += 1
                    }
                    cur_seq_idx += 1;
                }
                // assert that we extract the same number of modifications as we have distances
                assert_eq!(
                    cur_mod_idx,
                    mod_dists.len(),
                    "{:?} {}",
                    String::from_utf8_lossy(record.qname()),
                    record.is_reverse()
                );

                // check for the probability of modification.
                let num_mods_cur_end = num_mods_seen + unfiltered_modified_positions.len();
                let unfiltered_modified_probabilities = if num_mods_cur_end > ml_tag.len() {
                    let needed_num_of_zeros = num_mods_cur_end - ml_tag.len();
                    let mut to_add = vec![0; needed_num_of_zeros];
                    let mut has = ml_tag[num_mods_seen..ml_tag.len()].to_vec();
                    has.append(&mut to_add);
                    log::warn!(
                        "ML tag is too short for the number of modifications found in the MM tag. Assuming an ML value of 0 after the first {num_mods_cur_end} modifications."
                    );
                    has
                } else {
                    ml_tag[num_mods_seen..num_mods_cur_end].to_vec()
                };
                num_mods_seen = num_mods_cur_end;

                // must be true for filtering, and at this point
                assert_eq!(
                    unfiltered_modified_positions.len(),
                    unfiltered_modified_probabilities.len()
                );

                // Filter mods based on probabilities
                let (modified_probabilities, modified_positions): (Vec<u8>, Vec<i64>) =
                    unfiltered_modified_probabilities
                        .iter()
                        .zip(unfiltered_modified_positions.iter())
                        .filter(|(&ml, &_mm)| ml >= min_ml_score)
                        .unzip();

                // don't add empty basemods
                if modified_positions.is_empty() {
                    continue;
                }
                // add to a struct
                let mods = BaseMod::new(
                    record,
                    mod_base,
                    mod_strand.chars().next().unwrap(),
                    modification_type.chars().next().unwrap(),
                    modified_positions,
                    modified_probabilities,
                );
                rtn.push(mods);
            }
        } else {
            log::trace!("No MM tag found");
        }

        if ml_tag.len() != num_mods_seen {
            log::warn!(
                "ML tag ({}) different number than MM tag ({}).",
                ml_tag.len(),
                num_mods_seen
            );
        }
        // needed so I can compare methods
        rtn.sort();
        BaseMods { base_mods: rtn }
    }

    pub fn hashmap_to_basemods(
        map: HashMap<(i32, i32, i32), Vec<(i64, u8)>>,
        record: &bam::Record,
    ) -> BaseMods {
        let mut rtn = vec![];
        for (mod_info, mods) in map {
            let mod_base = mod_info.0 as u8;
            let mod_type = mod_info.1 as u8 as char;
            let mod_strand = if mod_info.2 == 0 { '+' } else { '-' };
            let (mut positions, mut qualities): (Vec<i64>, Vec<u8>) = mods.into_iter().unzip();
            if record.is_reverse() {
                let length = record.seq_len() as i64;
                positions = positions
                    .into_iter()
                    .rev()
                    .map(|p| length - p - 1)
                    .collect();
                qualities.reverse();
            }
            let mods = BaseMod::new(record, mod_base, mod_strand, mod_type, positions, qualities);
            rtn.push(mods);
        }
        // needed so I can compare methods
        rtn.sort();
        BaseMods { base_mods: rtn }
    }

    /// remove m6a base mods from the struct
    pub fn drop_m6a(&mut self) {
        self.base_mods.retain(|bm| !bm.is_m6a());
    }

    /// remove cpg/5mc base mods from the struct
    pub fn drop_cpg(&mut self) {
        self.base_mods.retain(|bm| !bm.is_cpg());
    }

    /// drop the forward stand of basemod calls
    pub fn drop_forward(&mut self) {
        self.base_mods.retain(|bm| bm.strand == '-');
    }

    /// drop the reverse strand of basemod calls   
    pub fn drop_reverse(&mut self) {
        self.base_mods.retain(|bm| bm.strand == '+');
    }

    /// drop m6A modifications with a qual less than the min_ml_score
    pub fn filter_m6a(&mut self, min_ml_score: u8) {
        self.base_mods
            .iter_mut()
            .filter(|bm| bm.is_m6a())
            .for_each(|bm| bm.ranges.filter_by_qual(min_ml_score));
    }

    /// drop 5mC modifications with a qual less than the min_ml_score
    pub fn filter_5mc(&mut self, min_ml_score: u8) {
        self.base_mods
            .iter_mut()
            .filter(|bm| bm.is_cpg())
            .for_each(|bm| bm.ranges.filter_by_qual(min_ml_score));
    }

    /// filter the basemods at the read ends
    pub fn filter_at_read_ends(&mut self, n_strip: i64) {
        if n_strip <= 0 {
            return;
        }
        self.base_mods
            .iter_mut()
            .for_each(|bm| bm.filter_at_read_ends(n_strip));
    }

    /// combine the forward and reverse m6a data
    pub fn m6a(&self) -> Ranges {
        let ranges = self
            .base_mods
            .iter()
            .filter(|x| x.is_m6a())
            .map(|x| &x.ranges)
            .collect();
        Ranges::merge_ranges(ranges)
    }

    /// combine the forward and reverse cpd/5mc data
    pub fn cpg(&self) -> Ranges {
        let ranges = self
            .base_mods
            .iter()
            .filter(|x| x.is_cpg())
            .map(|x| &x.ranges)
            .collect();
        Ranges::merge_ranges(ranges)
    }

    /// Example MM tag: MM:Z:C+m,11,6,10;A+a,0,0,0;
    /// Example ML tag: ML:B:C,157,30,2,164,118,255
    pub fn add_mm_and_ml_tags(&self, record: &mut bam::Record) {
        // init the mm and ml tag to be populated
        let mut ml_tag: Vec<u8> = vec![];
        let mut mm_tag = "".to_string();
        // need the original sequence for distances between bases.
        let mut seq = record.seq().as_bytes();
        if record.is_reverse() {
            seq = revcomp(seq);
        }
        // add to the ml and mm tag.
        for basemod in self.base_mods.iter() {
            // adding quality values (ML)
            ml_tag.extend(basemod.ranges.get_forward_quals());
            // get MM tag values
            let mut cur_mm = vec![];
            let positions = basemod.ranges.forward_starts();
            let mut last_pos = 0;
            for pos in positions {
                let u_pos = pos as usize;
                let mut in_between = 0;
                if last_pos < u_pos {
                    for base in seq[last_pos..u_pos].iter() {
                        if *base == basemod.modified_base {
                            in_between += 1;
                        }
                    }
                }
                last_pos = u_pos + 1;
                cur_mm.push(in_between);
            }
            // Add to the MM string
            mm_tag.push(basemod.modified_base as char);
            mm_tag.push(basemod.strand);
            mm_tag.push(basemod.modification_type);
            for diff in cur_mm {
                mm_tag.push_str(&format!(",{diff}"));
            }
            mm_tag.push(';')
            // next basemod
        }
        log::trace!(
            "{}\n{}\n{}\n",
            record.is_reverse(),
            mm_tag,
            String::from_utf8_lossy(&seq)
        );
        // clear out the old base mods
        record.remove_aux(b"MM").unwrap_or(());
        record.remove_aux(b"ML").unwrap_or(());
        // Add MM
        let aux_integer_field = Aux::String(&mm_tag);
        record.push_aux(b"MM", aux_integer_field).unwrap();
        // Add ML
        let aux_array: AuxArray<u8> = (&ml_tag).into();
        let aux_array_field = Aux::ArrayU8(aux_array);
        record.push_aux(b"ML", aux_array_field).unwrap();
    }
}

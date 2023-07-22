use bamlift::*;
use bio::alphabets::dna::revcomp;
use bio_io::*;
use itertools::{izip, multiunzip};
use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{
    bam,
    bam::record::{Aux, AuxArray},
};
use std::collections::HashMap;

use std::convert::TryFrom;

#[derive(Eq, PartialEq, Debug, PartialOrd, Ord)]
pub struct BaseMod {
    pub modified_base: u8,
    pub strand: char,
    pub modification_type: char,
    record_is_reverse: bool,
    modified_bases: Vec<i64>,
    modified_bases_forward: Vec<i64>,
    modified_probabilities: Vec<u8>,
    reference_positions: Vec<i64>,
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
        let modified_bases = positions_on_complimented_sequence(record, &modified_bases_forward);

        // rev qualities if needed
        let modified_probabilities = if record.is_reverse() {
            modified_probabilities_forward.into_iter().rev().collect()
        } else {
            modified_probabilities_forward
        };

        let record_is_reverse = record.is_reverse();

        // get the reference positions
        let reference_positions = lift_reference_positions_exact(record, &modified_bases);
        Self {
            modified_base,
            strand,
            modification_type,
            record_is_reverse,
            modified_bases,
            modified_bases_forward,
            modified_probabilities,
            reference_positions,
        }
    }

    pub fn get_reference_positions(&self) -> Vec<i64> {
        self.reference_positions.clone()
    }

    pub fn get_modified_bases(&self) -> Vec<i64> {
        self.modified_bases.clone()
    }

    pub fn get_modified_bases_forward(&self) -> Vec<i64> {
        self.modified_bases_forward.clone()
    }

    pub fn get_modified_probabilities(&self) -> Vec<u8> {
        self.modified_probabilities.clone()
    }

    pub fn get_modified_probabilities_forward(&self) -> Vec<u8> {
        if self.record_is_reverse {
            self.modified_probabilities.iter().rev().cloned().collect()
        } else {
            self.modified_probabilities.clone()
        }
    }

    pub fn is_m6a(&self) -> bool {
        self.modification_type == 'a'
    }

    pub fn is_cpg(&self) -> bool {
        self.modification_type == 'm'
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct BaseMods {
    pub base_mods: Vec<BaseMod>,
}

impl BaseMods {
    pub fn new(record: &bam::Record, min_ml_score: u8) -> BaseMods {
        // my basemod parser is ~25% faster than rust_htslib's
        // let new = BaseMods::rust_htslib_mm_ml_parser(record, min_ml_score);
        //assert_eq!(new, old);
        BaseMods::my_mm_ml_parser(record, min_ml_score)
    }

    pub fn my_mm_ml_parser(record: &bam::Record, min_ml_score: u8) -> BaseMods {
        // regex for matching the MM tag
        lazy_static! {
            static ref MM_RE: Regex =
                Regex::new(r"((([ACGTUN])([-+])([a-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
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
                    revcomp(record.seq().as_bytes())
                } else {
                    record.seq().as_bytes()
                };

                // find real positions in the forward sequence
                let mut cur_mod_idx = 0;
                let mut cur_seq_idx = 0;
                let mut dist_from_last_mod_base = 0;
                let mut unfiltered_modified_positions: Vec<i64> = vec![0; mod_dists.len()];
                while cur_seq_idx < forward_bases.len() && cur_mod_idx < mod_dists.len() {
                    let cur_base = forward_bases[cur_seq_idx];
                    if cur_base == mod_base && dist_from_last_mod_base == mod_dists[cur_mod_idx] {
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
                assert_eq!(cur_mod_idx, mod_dists.len());

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
            log::debug!("No MM tag found");
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

    /// this is actually way slower than my parser for some reason...
    /// so I am not going to use it
    pub fn rust_htslib_mm_ml_parser(record: &bam::Record, min_ml_score: u8) -> BaseMods {
        let mut base_mods = HashMap::new();

        match record.basemods_iter() {
            Ok(mods) => {
                for bp in mods {
                    let (pos, m) = bp.unwrap();
                    if min_ml_score > m.qual as u8 {
                        continue;
                    }
                    let z = (m.canonical_base, m.modified_base, m.strand);
                    let cur_mod = base_mods.entry(z).or_insert(vec![]);
                    cur_mod.push((pos as i64, m.qual as u8));
                    continue;
                }
            }
            _ => {
                return {
                    log::warn!("Rust hts-lib failed to parse basemods, trying custom parser.");
                    BaseMods::my_mm_ml_parser(record, min_ml_score)
                }
            }
        }
        BaseMods::hashmap_to_basemods(base_mods, record)
    }

    /// remove m6a base mods from the struct
    pub fn drop_m6a(&mut self) {
        self.base_mods.retain(|bm| !bm.is_m6a());
    }

    /// remove cpg/5mc base mods from the struct
    pub fn drop_cpg(&mut self) {
        self.base_mods.retain(|bm| !bm.is_cpg());
    }

    pub fn m6a_positions(&self, reference: bool) -> Vec<i64> {
        let m6a: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_m6a()).collect();
        // skip if no m6a
        if m6a.is_empty() {
            return vec![];
        }
        // if we only have one of the two mod types
        if m6a.len() == 1 {
            if reference {
                return m6a[0].get_reference_positions();
            } else {
                return m6a[0].get_modified_bases();
            }
        }
        // get positions of m6a if we have both A+a and T-a
        if reference {
            merge_two_lists(
                &m6a[0].get_reference_positions(),
                &m6a[1].get_reference_positions(),
            )
        } else {
            merge_two_lists(&m6a[0].get_modified_bases(), &m6a[1].get_modified_bases())
        }
    }

    pub fn m6a_full_probabilities(&self, record: &bam::Record) -> Vec<(i64, f32)> {
        let mp = get_f32_tag(record, b"mp");
        let m6a: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_m6a()).collect();
        // skip if no mp or m6a
        if m6a.is_empty() || mp.is_empty() {
            return vec![];
        }
        let m6a: Vec<i64> = m6a.iter().flat_map(|x| x.get_modified_bases()).collect();
        // skip if not equal
        if m6a.len() != mp.len() {
            log::warn!(
                "In {} m6A mods ({}) not equal to number of predictions ({}), returning nothing for this read.",
                String::from_utf8_lossy(record.qname()),
                m6a.len(),
                mp.len()
            );
            return vec![];
        }
        m6a.into_iter().zip(mp.into_iter()).collect()
    }

    fn helper_get_m6a(&self, forward: bool) -> (Vec<i64>, Vec<i64>, Vec<u8>) {
        let m6a: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_m6a()).collect();
        // skip if no m6a
        if m6a.is_empty() {
            return (vec![], vec![], vec![]);
        }

        let m6a_pos: Vec<i64> = if forward {
            m6a.iter()
                .flat_map(|x| x.get_modified_bases_forward())
                .collect()
        } else {
            m6a.iter().flat_map(|x| x.get_modified_bases()).collect()
        };

        let m6a_ref: Vec<i64> = m6a
            .iter()
            .flat_map(|x| x.get_reference_positions())
            .collect();
        let m6a_qual: Vec<u8> = m6a
            .iter()
            .flat_map(|x| x.get_modified_probabilities())
            .collect();

        assert_eq!(m6a_pos.len(), m6a_ref.len());
        assert_eq!(m6a_ref.len(), m6a_qual.len());
        let mut z: Vec<(i64, i64, u8)> = izip!(m6a_pos, m6a_ref, m6a_qual).collect();
        z.sort_by_key(|(p, _r, _q)| *p);
        multiunzip(z)
    }

    pub fn m6a(&self) -> (Vec<i64>, Vec<i64>, Vec<u8>) {
        self.helper_get_m6a(false)
    }

    pub fn forward_m6a(&self) -> (Vec<i64>, Vec<i64>, Vec<u8>) {
        self.helper_get_m6a(true)
    }

    pub fn cpg_positions(&self, reference: bool) -> Vec<i64> {
        let cpg: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_cpg()).collect();
        // skip if no cpg
        if cpg.is_empty() {
            return vec![];
        }
        // get positions of cpg
        if reference {
            cpg[0].get_reference_positions()
        } else {
            cpg[0].get_modified_bases()
        }
    }

    fn helper_get_cpg(&self, forward: bool) -> (Vec<i64>, Vec<i64>, Vec<u8>) {
        let cpg: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_cpg()).collect();
        // skip if no m6a
        if cpg.is_empty() {
            return (vec![], vec![], vec![]);
        }

        let cpg_pos: Vec<i64> = if forward {
            cpg.iter()
                .flat_map(|x| x.get_modified_bases_forward())
                .collect()
        } else {
            cpg.iter().flat_map(|x| x.get_modified_bases()).collect()
        };

        let cpg_ref: Vec<i64> = cpg
            .iter()
            .flat_map(|x| x.get_reference_positions())
            .collect();
        let cpg_qual: Vec<u8> = cpg
            .iter()
            .flat_map(|x| x.get_modified_probabilities())
            .collect();

        assert_eq!(cpg_pos.len(), cpg_ref.len());
        assert_eq!(cpg_ref.len(), cpg_qual.len());
        let mut z: Vec<(i64, i64, u8)> = izip!(cpg_pos, cpg_ref, cpg_qual).collect();
        z.sort_by_key(|(p, _r, _q)| *p);
        multiunzip(z)
    }

    pub fn cpg(&self) -> (Vec<i64>, Vec<i64>, Vec<u8>) {
        self.helper_get_cpg(false)
    }

    pub fn forward_cpg(&self) -> (Vec<i64>, Vec<i64>, Vec<u8>) {
        self.helper_get_cpg(true)
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
            ml_tag.extend(basemod.get_modified_probabilities_forward());
            // get MM tag values
            let mut cur_mm = vec![];
            let positions = basemod.get_modified_bases_forward();
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
                mm_tag.push_str(&format!(",{}", diff));
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

#[cfg(test)]
mod tests {
    use super::*;
    use env_logger::{Builder, Target};
    use log;
    use rust_htslib::{bam, bam::Read};

    #[test]
    /// checks that we can read base mods into BaseMods and the write the BaseMods to a new bam
    /// record, extract them again and they remain the same.
    fn test_mods_do_not_change() {
        Builder::new()
            .target(Target::Stderr)
            .filter(None, log::LevelFilter::Debug)
            .init();
        let mut bam = bam::Reader::from_path(&"tests/data/all.bam").unwrap();
        for rec in bam.records() {
            let mut rec = rec.unwrap();
            let mods = BaseMods::new(&rec, 0);
            mods.add_mm_and_ml_tags(&mut rec);
            let mods_2 = BaseMods::new(&rec, 0);
            assert_eq!(mods, mods_2);
        }
    }
}

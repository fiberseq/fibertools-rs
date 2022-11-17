use super::bamlift::*;
use bio::alphabets::dna::revcomp;
//use itertools::{izip, multiunzip};
use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{bam, bam::record::Aux};
use std::convert::TryFrom;

pub struct BaseMod {
    pub modified_base: u8,
    pub strand: char,
    pub modification_type: char,
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
        modified_probabilities: Vec<u8>,
    ) -> Self {
        let modified_bases = positions_on_complimented_sequence(record, &modified_bases_forward);
        // get the reference positions
        let reference_positions = get_exact_reference_positions(record, &modified_bases);
        Self {
            modified_base,
            strand,
            modification_type,
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
        if self.strand == '-' {
            self.modified_probabilities
                .clone()
                .into_iter()
                .rev()
                .collect()
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

pub struct BaseMods {
    pub base_mods: Vec<BaseMod>,
}

impl BaseMods {
    pub fn new(record: &bam::Record, min_ml_score: u8) -> BaseMods {
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

        if ml_tag.len() > num_mods_seen {
            log::warn!("ML tag has more entries than # of modifications in the MM tag.");
        }

        BaseMods { base_mods: rtn }
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

    pub fn m6a(&self) -> (Vec<i64>, Vec<i64>, Vec<u8>) {
        let m6a: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_m6a()).collect();
        // skip if no m6a
        if m6a.is_empty() {
            return (vec![], vec![], vec![]);
        }
        let m6a_pos: Vec<i64> = m6a.iter().flat_map(|x| x.get_modified_bases()).collect();
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
        //let mut _z: Vec<(i64, i64, u8)> = izip!(m6a_pos, m6a_ref, m6a_qual).collect();
        //_z.sort_by_key(|(p, _r, _q)| *p);
        //let _a = multiunzip(_z);
        (m6a_pos, m6a_ref, m6a_qual)
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

    /* TODO
    pub fn make_mm_and_ml_tags(&self) -> Vec<u8> {
        let mut mm = "".to_string();
        let mut ml = vec![];
        for basemod in &self.base_mods {
            for
        }
        ml
    }
    */
}

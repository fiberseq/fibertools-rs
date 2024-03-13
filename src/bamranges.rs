use super::bamlift::*;
use itertools::{izip, multiunzip};
use rust_htslib::bam;

#[derive(Debug, Clone, PartialEq, Eq, Ord, PartialOrd)]
pub struct Ranges {
    pub starts: Vec<Option<i64>>,
    pub ends: Vec<Option<i64>>,
    pub lengths: Vec<Option<i64>>,
    pub qual: Vec<u8>,
    pub reference_starts: Vec<Option<i64>>,
    pub reference_ends: Vec<Option<i64>>,
    pub reference_lengths: Vec<Option<i64>>,
    pub seq_len: i64,
    reverse: bool,
}

impl Ranges {
    /// starts and ends are [) intervals.
    pub fn new(
        record: &bam::Record,
        mut forward_starts: Vec<i64>,
        forward_ends: Option<Vec<i64>>,
        mut lengths: Option<Vec<i64>>,
    ) -> Self {
        let mut single_bp_liftover = false;
        // assume ends == starts if not provided
        if forward_ends.is_none() && lengths.is_none() {
            lengths = Some(vec![1; forward_starts.len()]);
            single_bp_liftover = true;
        }

        // use ends or calculate them
        let mut forward_ends_inclusive: Vec<i64> = match forward_ends {
            Some(x) => x.into_iter().map(|x| x + 1).collect(),
            None => forward_starts
                .iter()
                .zip(lengths.unwrap().iter())
                .map(|(&x, &y)| x + y - 1)
                .collect(),
        };

        // bam features for finding aligned positions
        let is_reverse = record.is_reverse();
        let seq_len = i64::try_from(record.seq_len()).unwrap();

        // get positions and lengths in reference orientation
        Self::positions_on_aligned_sequence(&mut forward_starts, is_reverse, seq_len);
        Self::positions_on_aligned_sequence(&mut forward_ends_inclusive, is_reverse, seq_len);
        let mut starts = forward_starts;
        let mut ends = forward_ends_inclusive;

        // swaps starts and ends if we reverse complemented
        if record.is_reverse() {
            std::mem::swap(&mut starts, &mut ends);
        }

        // swap back to non-inclusive ends
        ends = ends.into_iter().map(|x| x + 1).collect();

        // get lengths
        let lengths = starts
            .iter()
            .zip(ends.iter())
            .map(|(&x, &y)| Some(y - x))
            .collect::<Vec<_>>();

        let (reference_starts, reference_ends, reference_lengths) = if single_bp_liftover {
            lift_query_range_exact(record, &starts, &starts)
        } else {
            lift_query_range(record, &starts, &ends)
        };

        // return object
        Ranges {
            starts: starts.into_iter().map(Some).collect(),
            ends: ends.into_iter().map(Some).collect(),
            lengths,
            qual: vec![0; reference_starts.len()],
            reference_starts,
            reference_ends,
            reference_lengths,
            seq_len,
            reverse: is_reverse,
        }
    }

    pub fn set_qual(&mut self, qual: Vec<u8>) {
        assert_eq!(qual.len(), self.starts.len());
        self.qual = qual;
        if self.reverse {
            self.qual.reverse();
        }
    }

    /// get positions on the complimented sequence in the cigar record
    fn positions_on_aligned_sequence(input_positions: &mut [i64], is_reverse: bool, seq_len: i64) {
        if !is_reverse {
            return;
        }
        //need to correct for going from [) to (] if we are part of a range
        for p in input_positions.iter_mut() {
            *p = seq_len - *p - 1;
        }
        input_positions.reverse();
    }

    /// get the molecular coordinates of the ranges, taking into account
    /// the alignment orientation
    pub fn get_molecular(&self) -> Vec<Option<(i64, i64, i64)>> {
        self.starts
            .iter()
            .zip(self.ends.iter())
            .zip(self.lengths.iter())
            .map(|((s, e), l)| {
                if let (Some(s), Some(e), Some(l)) = (s, e, l) {
                    Some((*s, *e, *l))
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn get_starts(&self) -> Vec<i64> {
        self.starts.clone().into_iter().flatten().collect()
    }

    pub fn get_ends(&self) -> Vec<i64> {
        self.ends.clone().into_iter().flatten().collect()
    }

    pub fn get_forward_starts(&self) -> Vec<i64> {
        let mut z = self.get_starts();
        Self::positions_on_aligned_sequence(&mut z, self.reverse, self.seq_len);
        z
    }

    pub fn get_forward_quals(&self) -> Vec<u8> {
        let mut forward = self.qual.clone();
        if self.reverse {
            forward.reverse();
        }
        forward
    }

    pub fn to_strings(&self, reference: bool, skip_none: bool) -> Vec<String> {
        let (s, e, l, q) = if reference {
            (
                &self.reference_starts,
                &self.reference_ends,
                &self.reference_lengths,
                &self.qual,
            )
        } else {
            (&self.starts, &self.ends, &self.lengths, &self.qual)
        };

        let s = super::join_by_str_option_can_skip(s, ",", skip_none);
        let e = super::join_by_str_option_can_skip(e, ",", skip_none);
        let l = super::join_by_str_option_can_skip(l, ",", skip_none);
        if reference {
            vec![s, e, l]
        } else {
            let q = super::join_by_str(q, ",");
            vec![s, e, l, q]
        }
    }

    /// get the reference coordinates of the ranges, taking into account
    /// the alignment orientation
    pub fn get_reference(&self) -> Vec<Option<(i64, i64, i64)>> {
        self.reference_starts
            .iter()
            .zip(self.reference_ends.iter())
            .zip(self.reference_lengths.iter())
            .map(|((s, e), l)| {
                if let (Some(s), Some(e), Some(l)) = (s, e, l) {
                    Some((*s, *e, *l))
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn merge_ranges(multiple_ranges: Vec<&Self>) -> Self {
        if multiple_ranges.is_empty() {
            return Self {
                starts: vec![],
                ends: vec![],
                lengths: vec![],
                qual: vec![],
                reference_starts: vec![],
                reference_ends: vec![],
                reference_lengths: vec![],
                seq_len: 0,
                reverse: false,
            };
        }
        // check properties that must be the same
        let reverse = multiple_ranges[0].reverse;
        let seq_len = multiple_ranges[0].seq_len;
        for r in multiple_ranges.iter() {
            assert_eq!(r.reverse, reverse);
            assert_eq!(r.seq_len, seq_len);
        }
        // get the other properties
        let starts = multiple_ranges.iter().flat_map(|r| r.starts.clone());
        let ends = multiple_ranges.iter().flat_map(|r| r.ends.clone());
        let lengths = multiple_ranges.iter().flat_map(|r| r.lengths.clone());
        let qual = multiple_ranges.iter().flat_map(|r| r.qual.clone());
        let reference_starts = multiple_ranges
            .iter()
            .flat_map(|r| r.reference_starts.clone());
        let reference_ends = multiple_ranges
            .iter()
            .flat_map(|r| r.reference_ends.clone());
        let reference_lengths = multiple_ranges
            .iter()
            .flat_map(|r| r.reference_lengths.clone());

        #[allow(clippy::type_complexity)]
        let mut combo: Vec<(
            Option<i64>,
            Option<i64>,
            Option<i64>,
            u8,
            Option<i64>,
            Option<i64>,
            Option<i64>,
        )> = izip!(
            starts,
            ends,
            lengths,
            qual,
            reference_starts,
            reference_ends,
            reference_lengths
        )
        .collect();
        // sort by start position
        combo.sort_by_key(|(s, _e, _l, _q, _r_s, _r_e, _r_l)| *s);
        // unzip
        let (starts, ends, lengths, qual, reference_starts, reference_ends, reference_lengths) =
            multiunzip(combo);

        Self {
            starts,
            ends,
            lengths,
            qual,
            reference_starts,
            reference_ends,
            reference_lengths,
            seq_len,
            reverse,
        }
    }
}

impl<'a> IntoIterator for &'a Ranges {
    type Item = (i64, i64, i64, u8, Option<(i64, i64, i64)>);
    type IntoIter = RangesIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        RangesIterator {
            row: self,
            index: 0,
        }
    }
}

pub struct RangesIterator<'a> {
    row: &'a Ranges,
    index: usize,
}

impl<'a> Iterator for RangesIterator<'a> {
    type Item = (i64, i64, i64, u8, Option<(i64, i64, i64)>);
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.row.starts.len() {
            return None;
        }
        let start = self.row.starts[self.index]?;
        let end = self.row.ends[self.index]?;
        let length = self.row.lengths[self.index]?;
        let qual = self.row.qual[self.index];
        let reference_start = self.row.reference_starts[self.index];
        let reference_end = self.row.reference_ends[self.index];
        let reference_length = self.row.reference_lengths[self.index];
        let reference = match (reference_start, reference_end, reference_length) {
            (Some(s), Some(e), Some(l)) => Some((s, e, l)),
            _ => None,
        };
        self.index += 1;
        Some((start, end, length, qual, reference))
    }
}

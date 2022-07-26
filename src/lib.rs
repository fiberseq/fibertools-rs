/// Command line interface for fibertools-rs.
pub mod cli;

pub mod extract;

use rust_htslib::{bam, bam::ext::BamRecordExtensions};
use std::fmt::Display;
/// search a sorted array for insertions positions of another sorted array
/// returned index i satisfies
/// left
/// a[i-1] < v <= a[i]
/// right
/// a[i-1] <= v < a[i]
/// https://numpy.org/doc/stable/reference/generated/numpy.searchsorted.html
/// ```
/// use fibertools_rs::*;
/// let a = vec![1, 2, 3, 5, 6, 7, 8, 9, 10];
/// let v = vec![0, 1, 3, 4, 11, 11];
/// let indexes = search_sorted(&a, &v);
/// assert_eq!(indexes, vec![0, 0, 2, 3, 9, 9]);
/// ```
pub fn search_sorted<T>(a: &Vec<T>, v: &Vec<T>) -> Vec<usize>
where
    T: Ord,
    T: Display,
{
    let mut indexes = Vec::with_capacity(v.len());
    let mut a_idx = 0;
    for cur_v in v {
        while a_idx < a.len() {
            // check starting condition
            if a_idx == 0 && *cur_v <= a[a_idx] {
                indexes.push(0);
                break;
            } else if a_idx == 0 {
                a_idx += 1;
            }
            // end condition
            if a_idx == a.len() - 1 && *cur_v > a[a_idx] {
                indexes.push(a_idx + 1);
                break;
            }
            // middle of the array
            else if (a[a_idx - 1] < *cur_v) && (*cur_v <= a[a_idx]) {
                indexes.push(a_idx);
                break;
            }
            a_idx += 1;
        }
        //println!("{} {}\t{} {}", a_idx, v_idx, a[a_idx], v[v_idx]);
    }
    indexes
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
///     let seq_len = i64::try_from(record.seq_len()).unwrap();
///     let positions: Vec<i64> = (0..seq_len).collect();
///     liftover_closest(&record, &positions);    
/// }
///```
pub fn liftover_closest(record: &bam::Record, positions: &Vec<i64>) -> Vec<i64> {
    // find the shared positions in the reference
    let mut ref_positions = vec![];
    let (_q_pos, r_pos): (Vec<i64>, Vec<i64>) = record
        .aligned_pairs()
        .map(|[q_pos, r_pos]| (q_pos, r_pos))
        .unzip();
    let ref_idxs = search_sorted(&r_pos, positions);
    for mut idx in ref_idxs {
        // if we map past the end of the reference take the last reference position
        if idx == r_pos.len() {
            idx -= 1;
        }
        ref_positions.push(r_pos[idx]);
    }
    ref_positions
}

/// liftover positions using the cigar string
pub fn liftover_exact(record: &bam::Record, positions: &Vec<i64>) -> Vec<i64> {
    // find the shared positions in the reference
    let mut ref_positions = vec![];
    let mut cur_pos = 0;
    for [q_pos, r_pos] in record.aligned_pairs() {
        while cur_pos < positions.len() && positions[cur_pos] <= q_pos {
            if positions[cur_pos] == q_pos {
                log::trace!("Found position: q_pos:{}, r_pos:{}", q_pos, r_pos);
                ref_positions.push(r_pos);
            }
            cur_pos += 1;
        }
        if cur_pos == positions.len() {
            break;
        }
    }
    ref_positions
}

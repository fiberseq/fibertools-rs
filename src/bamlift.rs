use rust_htslib::{bam, bam::ext::BamRecordExtensions};
use std::fmt::Display;
/// Merge two lists into a sorted list
/// Normal sort is supposed to be very fast on two sorted lists
/// https://doc.rust-lang.org/std/vec/struct.Vec.html#current-implementation-6
pub fn merge_two_lists<T>(left: &[T], right: &[T]) -> Vec<T>
where
    T: Ord,
    T: Clone,
{
    let mut x: Vec<T> = left.iter().chain(right.iter()).cloned().collect();
    x.sort();
    x
}

/// get positions on the complimented sequence in the cigar record
pub fn positions_on_complimented_sequence(
    record: &bam::Record,
    input_positions: &[i64],
) -> Vec<i64> {
    // reverse positions if needed
    let positions: Vec<i64> = if record.is_reverse() {
        let seq_len = i64::try_from(record.seq_len()).unwrap();
        input_positions.iter().rev().map(|p| seq_len - p).collect()
    } else {
        input_positions.to_vec()
    };
    positions
}

/// search a sorted array for insertions positions of another sorted array
/// returned index i satisfies
/// left
/// a[i-1] < v <= a[i]
/// right
/// a[i-1] <= v < a[i]
/// https://numpy.org/doc/stable/reference/generated/numpy.searchsorted.html
/// ```
/// use fibertools_rs::bamlift::*;
/// let a = vec![1, 2, 3, 5, 6, 7, 8, 9, 10];
/// let v = vec![0, 1, 3, 4, 11, 11];
/// let indexes = search_sorted(&a, &v);
/// assert_eq!(indexes, vec![0, 0, 2, 3, 9, 9]);
/// ```
pub fn search_sorted<T>(a: &[T], v: &[T]) -> Vec<usize>
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
/// use fibertools_rs::bamlift::*;
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
pub fn liftover_closest(record: &bam::Record, query_positions: &[i64]) -> Vec<i64> {
    // aligned pairs
    let (q_pos, r_pos): (Vec<i64>, Vec<i64>) = record
        .aligned_pairs()
        .map(|[q_pos, r_pos]| (q_pos, r_pos))
        .unzip();
    // find the closest position within the q_pos matches
    let ref_idxs = search_sorted(&q_pos, query_positions);
    // find the shared positions in the reference
    let mut ref_positions = vec![];
    for mut idx in ref_idxs {
        // if we map past the end of the reference take the last reference position
        if idx == r_pos.len() {
            idx -= 1;
        }
        //log::trace!("Idx {}\tr_pos {}", idx, r_pos[idx]);
        ref_positions.push(r_pos[idx]);
    }
    ref_positions
}

/// liftover positions using the cigar string
pub fn liftover_exact(record: &bam::Record, positions: &[i64]) -> Vec<i64> {
    // find the shared positions in the reference
    let mut ref_positions = vec![];
    let mut cur_pos = 0;
    for [q_pos, r_pos] in record.aligned_pairs() {
        while cur_pos < positions.len() && positions[cur_pos] <= q_pos {
            if positions[cur_pos] == q_pos {
                //log::trace!("Found position: q_pos:{}, r_pos:{}", q_pos, r_pos);
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

pub fn get_closest_reference_positions(positions: &[i64], record: &bam::Record) -> Vec<i64> {
    //
    //log::trace!("seq {:?}", positions);
    let positions = positions_on_complimented_sequence(record, positions);
    log::trace!("seq comp {:?}", positions);
    // get the reference positions
    liftover_closest(record, &positions)
}

pub fn get_closest_reference_range(
    starts: &[i64],
    lengths: &[i64],
    record: &bam::Record,
) -> Vec<(i64, i64)> {
    let mol_ends: Vec<i64> = starts
        .iter()
        .zip(lengths.iter())
        .map(|(start, length)| start + length)
        .collect();
    //log::trace!("seq {:?}", starts);
    //log::trace!("seq len {:?}", mol_ends);
    let ref_starts = get_closest_reference_positions(starts, record);
    let ref_ends = get_closest_reference_positions(&mol_ends, record);
    assert_eq!(ref_starts.len(), ref_ends.len());
    //log::trace!("ref_starts: {:?}", ref_starts);
    //log::trace!("ref_ends: {:?}", ref_ends);
    ref_starts
        .iter()
        .zip(ref_ends.iter())
        .filter(|(&start, &end)| start <= end) // filter out zero length ranges, basically means there is no liftover
        .map(|(&start, &end)| (start, end - start))
        .collect()
}

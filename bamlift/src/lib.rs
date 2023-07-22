use itertools::multiunzip;
use rust_htslib::{bam, bam::ext::BamRecordExtensions};
use std::collections::HashMap;
use std::fmt::{Debug, Display};

/// Merge two lists into a sorted list
/// Normal sort is supposed to be very fast on two sorted lists
/// <https://doc.rust-lang.org/std/vec/struct.Vec.html#current-implementation-6>
pub fn merge_two_lists<T>(left: &[T], right: &[T]) -> Vec<T>
where
    T: Ord,
    T: Clone,
{
    let mut x: Vec<T> = left.iter().chain(right.iter()).cloned().collect();
    x.sort();
    x
}
/// Merge two lists based on a key
/// Normal sort is supposed to be very fast on two sorted lists
/// <https://doc.rust-lang.org/std/vec/struct.Vec.html#current-implementation-6>
/// ```
/// use bamlift::*;
/// let x = vec![1,3];
/// let x_q = vec!["a","b"];
/// let y = vec![2,4];
/// let y_q = vec!["c", "d"];
/// let z = merge_two_lists_with_qual(&x, &x_q, &y, &y_q);
/// assert_eq!(z, vec![(1,"a"), (2,"c"), (3,"b"), (4, "d")]);
/// ```
pub fn merge_two_lists_with_qual<T, U>(
    left: &[T],
    left_q: &[U],
    right: &[T],
    right_q: &[U],
) -> Vec<(T, U)>
where
    T: Ord,
    T: Clone,
    U: Clone,
{
    let l = left
        .iter()
        .zip(left_q.iter())
        .map(|(a, b)| (a.clone(), b.clone()));

    let r = right
        .iter()
        .zip(right_q.iter())
        .map(|(a, b)| (a.clone(), b.clone()));

    let mut x: Vec<(T, U)> = l.chain(r).collect();
    x.sort_by_key(|(a, _b)| a.clone());
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
        input_positions
            .iter()
            .rev()
            .map(|p| seq_len - p - 1)
            .collect()
    } else {
        input_positions.to_vec()
    };
    positions
}

/// get positions on the complimented sequence in the cigar record
pub fn positions_on_complimented_sequence_in_place(
    record: &bam::Record,
    input_positions: &mut Vec<i64>,
    part_of_range: bool,
) {
    if !record.is_reverse() {
        return;
    }
    let seq_len = i64::try_from(record.seq_len()).unwrap();
    // need to correct for going from [) to (] if we are part of a range
    let offset = if part_of_range { 0 } else { 1 };
    for p in input_positions.iter_mut() {
        *p = seq_len - *p - offset;
    }
    input_positions.reverse();
}

#[inline(always)]
pub fn is_sorted<T>(v: &[T]) -> bool
where
    T: Ord,
{
    v.windows(2).all(|w| w[0] <= w[1])
}

/// search a sorted array for insertions positions of another sorted array
/// returned index i satisfies
/// left
/// a\[i-1\] < v <= a\[i\]
/// right
/// a\[i-1\] <= v < a\[i\]
/// <https://numpy.org/doc/stable/reference/generated/numpy.searchsorted.html>
/// ```
/// use bamlift::*;
/// let a = vec![1, 2, 3, 5, 6, 7, 8, 9, 10];
/// let v = vec![0, 1, 3, 4, 11, 11];
/// let indexes = search_sorted(&a, &v);
/// assert_eq!(indexes, vec![0, 0, 2, 3, 9, 9]);
/// ```
pub fn search_sorted<T>(a: &[T], v: &[T]) -> Vec<usize>
where
    T: Ord,
    T: Display,
    [T]: Debug,
{
    if !is_sorted(v) {
        panic!("v is not sorted: {:?}", v);
    }

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
    }
    log::trace!("search_sorted: {:?}\n{:?}", v, indexes);
    indexes
}

/// this is a helper function for liftover_closest that should only be called from there
/// The exception for this is test cases, where it should be easier to test this function
/// directly.
fn liftover_closest(
    positions: &[i64],
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
) -> Vec<i64> {
    // skip empty
    if positions.is_empty() {
        return vec![];
    }
    if aligned_block_pairs.is_empty() {
        return positions.iter().map(|_x| -1).collect();
    }
    assert!(is_sorted(positions));

    // find the closest position for every position
    let mut starting_block = 0;
    let mut pos_mapping = HashMap::new();
    for cur_pos in positions {
        pos_mapping.insert(cur_pos, (-1, i64::MAX));
        let mut current_block = 0;
        for ([q_st, q_en], [r_st, r_en]) in &aligned_block_pairs[starting_block..] {
            // get the previous closest position
            let (best_r_pos, best_diff) = pos_mapping.get_mut(cur_pos).unwrap();
            // exact match found
            if cur_pos >= &q_st && cur_pos < &q_en {
                let dist_from_start = cur_pos - q_st;
                *best_diff = 0;
                *best_r_pos = r_st + dist_from_start;
                break;
            }
            // we are before the start of the block
            else if cur_pos < &q_st {
                let diff = (q_st - cur_pos).abs();
                if diff < *best_diff {
                    *best_diff = diff;
                    *best_r_pos = *r_st;
                }
            }
            // we are past the end of the block
            else if cur_pos >= &q_en {
                let diff = (q_en - cur_pos).abs();
                if diff < *best_diff {
                    *best_diff = diff;
                    *best_r_pos = *r_en;
                }
                // we don't need to return to previous blocks since the input is sorted
                starting_block = current_block;
            }
            current_block += 1;
        }
    }
    let mut rtn = vec![];
    for q_pos in positions {
        let (r_pos, _diff) = pos_mapping.get(q_pos).unwrap();
        rtn.push(*r_pos);
    }
    assert_eq!(rtn.len(), positions.len());
    rtn
}

/// find the closest reference positions for a list of query positions
pub fn lift_reference_positions(
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
    query_positions: &[i64],
) -> Vec<i64> {
    liftover_closest(query_positions, aligned_block_pairs)
}

/// find the closest query positions for a list of reference positions
pub fn lift_query_positions(
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
    reference_positions: &[i64],
) -> Vec<i64> {
    // if lifting to the query, we need to reverse the pairs
    let aligned_block_pairs = aligned_block_pairs.iter().map(|(q, r)| (*r, *q)).collect();
    liftover_closest(reference_positions, &aligned_block_pairs)
}

pub fn lift_range(
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
    starts: &[i64],
    ends: &[i64],
    lift_reference_to_query: bool,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    assert_eq!(starts.len(), ends.len());
    let (ref_starts, ref_ends) = if !lift_reference_to_query {
        (
            lift_reference_positions(aligned_block_pairs, starts),
            lift_reference_positions(aligned_block_pairs, ends),
        )
    } else {
        (
            lift_query_positions(aligned_block_pairs, starts),
            lift_query_positions(aligned_block_pairs, ends),
        )
    };
    assert_eq!(ref_starts.len(), ref_ends.len());
    let rtn = ref_starts
        .into_iter()
        .zip(ref_ends.into_iter())
        .map(|(start, end)| {
            if end <= start {
                (-1, -1, -1)
            } else {
                (start, end, end - start)
            }
        })
        .collect::<Vec<_>>();
    multiunzip(rtn)
}

/// Find the closest range but hopefully better
pub fn lift_query_range(
    record: &bam::Record,
    starts: &[i64],
    ends: &[i64],
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // get the aligned block pairs
    let aligned_block_pairs: Vec<([i64; 2], [i64; 2])> = record.aligned_block_pairs().collect();
    lift_range(&aligned_block_pairs, starts, ends, false)
}

/// liftover positions using the cigar string
pub fn liftover_exact(
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
    positions: &[i64],
    lift_reference_to_query: bool,
) -> Vec<i64> {
    // find the shared positions in the reference
    let mut return_positions = vec![];
    let mut cur_idx = 0;
    // ends are not inclusive, I checked.
    for ([q_st, q_en], [r_st, r_en]) in aligned_block_pairs {
        let (st, en) = if !lift_reference_to_query {
            (q_st, q_en)
        } else {
            (r_st, r_en)
        };
        // check bounds
        if cur_idx == positions.len() {
            break;
        }
        let mut cur_pos = positions[cur_idx];
        // need to go to the next block
        while cur_pos < *en {
            if cur_pos >= *st {
                let dist_from_start = cur_pos - st;
                let rtn_pos = if !lift_reference_to_query {
                    r_st + dist_from_start
                } else {
                    q_st + dist_from_start
                };
                return_positions.push(rtn_pos);
            } else {
                return_positions.push(-1);
            }
            // reset current position
            cur_idx += 1;
            if cur_idx == positions.len() {
                break;
            }
            cur_pos = positions[cur_idx];
        }
    }

    // add values for things that won't lift at the end
    while positions.len() > return_positions.len() {
        return_positions.push(-1);
    }
    assert_eq!(positions.len(), return_positions.len());
    return_positions
}

pub fn lift_reference_positions_exact(record: &bam::Record, query_positions: &[i64]) -> Vec<i64> {
    if record.is_unmapped() {
        query_positions.iter().map(|_x| -1).collect()
    } else {
        let aligned_block_pairs: Vec<([i64; 2], [i64; 2])> = record.aligned_block_pairs().collect();
        liftover_exact(&aligned_block_pairs, query_positions, false)
    }
}

pub fn lift_query_positions_exact(record: &bam::Record, reference_positions: &[i64]) -> Vec<i64> {
    if record.is_unmapped() {
        reference_positions.iter().map(|_x| -1).collect()
    } else {
        let aligned_block_pairs: Vec<([i64; 2], [i64; 2])> = record.aligned_block_pairs().collect();
        liftover_exact(&aligned_block_pairs, reference_positions, true)
    }
}

use anyhow::Result;
use itertools::multiunzip;
use rust_htslib::{bam, bam::ext::BamRecordExtensions};
use std::collections::HashMap;

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
/// use fibertools_rs::utils::bamlift::*;
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

#[inline(always)]
pub fn is_sorted<T>(v: &[T]) -> bool
where
    T: Ord,
{
    v.windows(2).all(|w| w[0] <= w[1])
}

//
// CLOSEST LIFTOVER FUNCTIONS
//

/// this is a helper function for liftover_closest that should only be called from there
/// The exception for this is test cases, where it should be easier to test this function
/// directly.
#[allow(clippy::all)]
fn liftover_closest(
    positions: &[i64],
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
) -> Result<Vec<Option<i64>>> {
    // skip empty
    if positions.is_empty() {
        return Ok(vec![]);
    }
    if aligned_block_pairs.is_empty() {
        return Ok(positions.iter().map(|_x| None).collect());
    }
    if !is_sorted(positions) {
        return Err(anyhow::anyhow!(
            "Positions must be sorted before calling liftover! Array length: {}, first 5: {:?}, last 5: {:?}",
            positions.len(),
            &positions[..std::cmp::min(5, positions.len())],
            &positions[positions.len().saturating_sub(5)..]
        ));
    }

    // find the closest position for every position
    let mut starting_block = 0;
    let ending_block = aligned_block_pairs.len();
    let mut pos_mapping = HashMap::new();
    for cur_pos in positions {
        pos_mapping.insert(cur_pos, (-1, i64::MAX));
        let mut current_block = 0;
        for block_index in starting_block..ending_block {
            // get the current alignment block
            let ([q_st, q_en], [r_st, r_en]) = &aligned_block_pairs[block_index];
            // get the previous closest position
            let (best_r_pos, best_diff) = pos_mapping.get_mut(cur_pos).unwrap();
            // exact match found (including exact end of block)
            if cur_pos >= q_st && cur_pos <= q_en {
                let dist_from_start = cur_pos - q_st;
                *best_diff = 0;
                *best_r_pos = r_st + dist_from_start;
                break;
            }
            // we are before the start of the block
            else if cur_pos < q_st {
                let diff = (q_st - cur_pos).abs();
                if diff < *best_diff {
                    *best_diff = diff;
                    *best_r_pos = *r_st;
                }
            }
            // we are past the end of the block
            else if cur_pos > q_en {
                let diff = (cur_pos - q_en).abs();
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
        let (r_pos, diff) = pos_mapping.get(q_pos).unwrap();
        if *r_pos == -1 && *diff == i64::MAX {
            rtn.push(None);
        } else {
            rtn.push(Some(*r_pos));
        }
    }
    assert_eq!(rtn.len(), positions.len());
    Ok(rtn)
}

/// find the closest reference positions for a list of query positions
pub fn lift_reference_positions(
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
    query_positions: &[i64],
) -> Result<Vec<Option<i64>>> {
    liftover_closest(query_positions, aligned_block_pairs)
}

/// find the closest query positions for a list of reference positions
pub fn lift_query_positions(
    aligned_block_pairs: &[([i64; 2], [i64; 2])],
    reference_positions: &[i64],
) -> Result<Vec<Option<i64>>> {
    // if lifting to the query, we need to reverse the pairs
    let aligned_block_pairs = aligned_block_pairs.iter().map(|(q, r)| (*r, *q)).collect();
    liftover_closest(reference_positions, &aligned_block_pairs)
}

/// Helper function to lift positions that may or may not be sorted
fn lift_positions_with_sort_handling(
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
    positions: &[i64],
    lift_reference_to_query: bool,
) -> Result<Vec<Option<i64>>> {
    if is_sorted(positions) {
        // Fast path: positions are sorted, lift normally
        if !lift_reference_to_query {
            lift_reference_positions(aligned_block_pairs, positions)
        } else {
            lift_query_positions(aligned_block_pairs, positions)
        }
    } else {
        // Slow path: positions are unsorted, need to sort positions independently
        let mut indices: Vec<usize> = (0..positions.len()).collect();
        indices.sort_by_key(|&i| positions[i]);

        let sorted_positions: Vec<i64> = indices.iter().map(|&i| positions[i]).collect();

        // Lift sorted positions
        let lifted_positions = if !lift_reference_to_query {
            lift_reference_positions(aligned_block_pairs, &sorted_positions)?
        } else {
            lift_query_positions(aligned_block_pairs, &sorted_positions)?
        };

        // Restore original order
        let mut original_positions = vec![None; positions.len()];
        for (sorted_idx, &original_idx) in indices.iter().enumerate() {
            original_positions[original_idx] = lifted_positions[sorted_idx];
        }

        Ok(original_positions)
    }
}

#[allow(clippy::type_complexity)]
fn lift_range(
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
    starts: &[i64],
    ends: &[i64],
    lift_reference_to_query: bool,
) -> Result<(Vec<Option<i64>>, Vec<Option<i64>>, Vec<Option<i64>>)> {
    assert_eq!(starts.len(), ends.len());

    // Lift starts and ends using the helper function
    let ref_starts =
        lift_positions_with_sort_handling(aligned_block_pairs, starts, lift_reference_to_query)?;
    let ref_ends =
        lift_positions_with_sort_handling(aligned_block_pairs, ends, lift_reference_to_query)?;

    // Common logic for processing lifted positions
    assert_eq!(ref_starts.len(), ref_ends.len());
    let rtn = ref_starts
        .into_iter()
        .zip(ref_ends)
        .map(|(start, end)| match (start, end) {
            (Some(start), Some(end)) => {
                if start == end {
                    (None, None, None)
                } else {
                    (Some(start), Some(end), Some(end - start))
                }
            }
            _ => (None, None, None),
        })
        .collect::<Vec<_>>();
    Ok(multiunzip(rtn))
}

/// Find the closest range but hopefully better
#[allow(clippy::type_complexity)]
pub fn lift_query_range(
    record: &bam::Record,
    starts: &[i64],
    ends: &[i64],
) -> Result<(Vec<Option<i64>>, Vec<Option<i64>>, Vec<Option<i64>>)> {
    // get the aligned block pairs
    let aligned_block_pairs: Vec<([i64; 2], [i64; 2])> = record.aligned_block_pairs().collect();
    lift_range(&aligned_block_pairs, starts, ends, false)
}

//
// EXACT LIFTOVER FUNCTIONS
//

/// liftover positions using the cigar string
fn liftover_exact(
    aligned_block_pairs: &Vec<([i64; 2], [i64; 2])>,
    positions: &[i64],
    lift_reference_to_query: bool,
) -> Result<Vec<Option<i64>>> {
    if !is_sorted(positions) {
        return Err(anyhow::anyhow!(
            "Positions must be sorted before calling liftover!"
        ));
    }

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
                return_positions.push(Some(rtn_pos));
            } else {
                return_positions.push(None);
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
        return_positions.push(None);
    }
    assert_eq!(positions.len(), return_positions.len());
    Ok(return_positions)
}

pub fn lift_reference_positions_exact(
    record: &bam::Record,
    query_positions: &[i64],
) -> Result<Vec<Option<i64>>> {
    if record.is_unmapped() {
        Ok(query_positions.iter().map(|_x| None).collect())
    } else {
        let aligned_block_pairs: Vec<([i64; 2], [i64; 2])> = record.aligned_block_pairs().collect();
        liftover_exact(&aligned_block_pairs, query_positions, false)
    }
}

pub fn lift_query_positions_exact(
    record: &bam::Record,
    reference_positions: &[i64],
) -> Result<Vec<Option<i64>>> {
    if record.is_unmapped() {
        Ok(reference_positions.iter().map(|_x| None).collect())
    } else {
        let aligned_block_pairs: Vec<([i64; 2], [i64; 2])> = record.aligned_block_pairs().collect();
        liftover_exact(&aligned_block_pairs, reference_positions, true)
    }
}

#[allow(clippy::type_complexity)]
fn lift_range_exact(
    record: &bam::Record,
    starts: &[i64],
    ends: &[i64],
    lift_reference_to_query: bool,
) -> Result<(Vec<Option<i64>>, Vec<Option<i64>>, Vec<Option<i64>>)> {
    assert_eq!(starts.len(), ends.len());
    let (ref_starts, ref_ends) = if !lift_reference_to_query {
        (
            lift_reference_positions_exact(record, starts)?,
            lift_reference_positions_exact(record, ends)?,
        )
    } else {
        (
            lift_query_positions_exact(record, starts)?,
            lift_query_positions_exact(record, ends)?,
        )
    };
    assert_eq!(ref_starts.len(), ref_ends.len());
    let rtn = ref_starts
        .into_iter()
        .zip(ref_ends)
        .map(|(start, end)| match (start, end) {
            (Some(start), Some(end)) => (Some(start), Some(end), Some(end - start)),
            _ => (None, None, None),
        })
        .collect::<Vec<_>>();
    Ok(multiunzip(rtn))
}

/// Find the exact range in the reference from a query range
#[allow(clippy::type_complexity)]
pub fn lift_query_range_exact(
    record: &bam::Record,
    starts: &[i64],
    ends: &[i64],
) -> Result<(Vec<Option<i64>>, Vec<Option<i64>>, Vec<Option<i64>>)> {
    lift_range_exact(record, starts, ends, false)
}

/// Find the exact range in the query from a reference range
#[allow(clippy::type_complexity)]
pub fn lift_reference_range_exact(
    record: &bam::Record,
    starts: &[i64],
    ends: &[i64],
) -> Result<(Vec<Option<i64>>, Vec<Option<i64>>, Vec<Option<i64>>)> {
    lift_range_exact(record, starts, ends, true)
}

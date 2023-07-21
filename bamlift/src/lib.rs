use rust_htslib::{bam, bam::ext::BamRecordExtensions};
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

fn liftover_closest(record: &bam::Record, positions: &[i64], get_reference: bool) -> Vec<i64> {
    // skip unmapped
    if record.is_unmapped() {
        return positions.iter().map(|_x| -1).collect();
    }

    // real work
    let (q_pos, r_pos): (Vec<i64>, Vec<i64>) = record
        .aligned_pairs()
        .map(|[q_pos, r_pos]| (q_pos, r_pos))
        .unzip();
    // if pulling from the query, we need to reverse the positions
    let (query_positions, return_positions) = if get_reference {
        (q_pos, r_pos)
    } else {
        (r_pos, q_pos)
    };

    // find the closest position within the q_pos matches
    let return_idxs = search_sorted(&query_positions, positions);
    // find the shared positions in the return positions
    let mut rtn = vec![];
    for mut idx in return_idxs {
        // if we map past the end of the reference take the last reference position
        if idx == return_positions.len() {
            idx -= 1;
        }
        //log::trace!("Idx {}\tr_pos {}", idx, r_pos[idx]);
        rtn.push(return_positions[idx]);
    }
    assert_eq!(rtn.len(), positions.len());
    rtn
}

///```
/// use rust_htslib::{bam, bam::Read};
/// use bamlift::*;
/// use log;
/// use env_logger::{Builder, Target};;
/// Builder::new().target(Target::Stderr).filter(None, log::LevelFilter::Debug).init();
/// let mut bam = bam::Reader::from_path(&"../tests/data/all.bam").unwrap();
/// for record in bam.records() {
///     let record = record.unwrap();
///     let seq_len = i64::try_from(record.seq_len()).unwrap();
///     let positions: Vec<i64> = (0..seq_len).collect();
///     closest_reference_positions(&record, &positions);
/// }
///```
pub fn closest_reference_positions(record: &bam::Record, query_positions: &[i64]) -> Vec<i64> {
    liftover_closest(record, query_positions, true)
}

pub fn closest_query_positions(record: &bam::Record, query_positions: &[i64]) -> Vec<i64> {
    liftover_closest(record, query_positions, false)
}

pub fn get_closest_reference_range(
    starts: &[i64],
    lengths: &[i64],
    record: &bam::Record,
) -> Vec<(i64, i64)> {
    assert_eq!(starts.len(), lengths.len());
    let ends: Vec<i64> = starts
        .iter()
        .zip(lengths.iter())
        .map(|(start, length)| start + length)
        .collect();

    let ref_starts = closest_reference_positions(record, starts);
    let ref_ends = closest_reference_positions(record, &ends);
    assert_eq!(ref_starts.len(), ref_ends.len());

    ref_starts
        .iter()
        .zip(ref_ends.iter())
        //.filter(|(&start, &end)| start < end) // filter out zero length ranges, basically means there is no liftover
        //.map(|(&start, &end)| (start, end - start))
        .map(|(&start, &end)| {
            if end <= start {
                (-1, -1)
            } else {
                (start, end - start)
            }
        })
        .collect()
}

/// liftover positions using the cigar string
fn _liftover_exact_old(record: &bam::Record, positions: &[i64], get_reference: bool) -> Vec<i64> {
    // find the shared positions in the reference
    let mut return_positions = vec![];
    let mut cur_pos = 0;
    for [q_pos, r_pos] in record.aligned_pairs() {
        // check whether we are searching for an exact position in
        // the reference or the query
        let val_to_match = if get_reference { q_pos } else { r_pos };
        // iterate over positions until we find the exact position or move past it
        while cur_pos < positions.len() && positions[cur_pos] <= val_to_match {
            if positions[cur_pos] == val_to_match {
                if get_reference {
                    return_positions.push(r_pos);
                } else {
                    return_positions.push(q_pos);
                }
            } else {
                return_positions.push(-1);
            }
            cur_pos += 1;
        }
        if cur_pos == positions.len() {
            break;
        }
    }
    // add values for things that won't lift at the end
    while positions.len() > return_positions.len() {
        return_positions.push(-1);
    }
    assert_eq!(positions.len(), return_positions.len());
    return_positions
}

/// liftover positions using the cigar string
fn liftover_exact_new(record: &bam::Record, positions: &[i64], get_reference: bool) -> Vec<i64> {
    // find the shared positions in the reference
    let mut return_positions = vec![];
    let mut cur_idx = 0;
    // ends are not inclusive, I checked.
    for ([q_st, q_en], [r_st, r_en]) in record.aligned_block_pairs() {
        let (st, en) = if get_reference {
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
        while cur_pos < en {
            if cur_pos >= st {
                let dist_from_start = cur_pos - st;
                let rtn_pos = if get_reference {
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

fn liftover_exact(record: &bam::Record, positions: &[i64], get_reference: bool) -> Vec<i64> {
    let new = liftover_exact_new(record, positions, get_reference);
    //let old = _liftover_exact_old(record, positions, get_reference);
    //assert_eq!(new, old);
    new
}

pub fn get_exact_reference_positions(record: &bam::Record, query_positions: &[i64]) -> Vec<i64> {
    if record.is_unmapped() {
        query_positions.iter().map(|_x| -1).collect()
    } else {
        liftover_exact(record, query_positions, true)
    }
}

pub fn get_exact_query_positions(record: &bam::Record, reference_positions: &[i64]) -> Vec<i64> {
    if record.is_unmapped() {
        reference_positions.iter().map(|_x| -1).collect()
    } else {
        liftover_exact(record, reference_positions, false)
    }
}

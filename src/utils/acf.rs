// Code browsed from:
// https://github.com/krfricke/arima
use anyhow::Result;
use num::Float;
use rayon::prelude::*;
use std::cmp;
use std::convert::From;
use std::ops::{Add, AddAssign, Div};

/// Calculate the auto-correlation function of a time series of length n.
///
/// # Arguments
///
/// * `&x` - Reference to input vector slice of length n.
/// * `max_lag` - Calculate ACF for this maximum lag. Defaults to n-1.
/// * `covariance` - If true, returns auto-covariances. If false, returns auto-correlations.
///
/// # Returns
///
/// * Output vector of length max_lag+1.
///
/// # Example
///
/// ```
/// use fibertools_rs::utils::acf;
/// let x = [1.0, 1.2, 1.4, 1.6];
/// let ac = acf::acf(&x, Some(2), false).unwrap();
/// assert!((ac[0] - 1.0) * (ac[0] - 1.0) < 1.0e-14);
/// assert!((ac[1] - 0.25) * (ac[1] - 0.25) < 1.0e-14);
/// assert!((ac[2] - (-0.3)) * (ac[2] - (-0.3)) < 1.0e-14);
/// ```
pub fn acf<T: Float + From<u32> + From<f64> + Copy + Add + AddAssign + Div>(
    x: &[T],
    max_lag: Option<usize>,
    covariance: bool,
) -> Result<Vec<T>> {
    let max_lag = match max_lag {
        // if upper bound for max_lag is n-1
        Some(max_lag) => cmp::min(max_lag, x.len() - 1),
        None => x.len() - 1,
    };
    if x.len() <= max_lag {
        return Err(anyhow::anyhow!(
            "acf-max-lag ({}) must be less than the number of m6A observations ({}).",
            max_lag,
            x.len()
        ));
    }

    let m = max_lag + 1;

    let len_x_usize = x.len();
    let len_x: T = From::from(len_x_usize as u32);
    let sum: T = From::from(0.0);

    let sum_x: T = x.iter().fold(sum, |sum, &xi| sum + xi);
    let mean_x: T = sum_x / len_x;

    let mut y: Vec<T> = vec![From::from(0.0); m];

    for t in 0..m {
        for i in 0..len_x_usize - t {
            let xi = x[i] - mean_x;
            let xi_t = x[i + t] - mean_x;
            y[t] += (xi * xi_t) / len_x;
        }
        // we need y[0] to calculate the correlations, so we set it to 1.0 at the end
        if !covariance && t > 0 {
            y[t] = y[t] / y[0];
        }
    }
    if !covariance {
        y[0] = From::from(1.0);
    }
    Ok(y)
}

/// Calculate the auto-correlation function of a time series of length n.
/// but this version is multithreaded over m using rayon.
pub fn acf_par<
    T: Float
        + From<u32>
        + From<f64>
        + Copy
        + Add
        + AddAssign
        + Div
        + Send
        + std::iter::Sum
        + std::marker::Sync,
>(
    x: &[T],
    max_lag: Option<usize>,
    covariance: bool,
) -> Result<Vec<T>> {
    let max_lag = match max_lag {
        // if upper bound for max_lag is n-1
        Some(max_lag) => cmp::min(max_lag, x.len() - 1),
        None => x.len() - 1,
    };
    if x.len() <= max_lag {
        return Err(anyhow::anyhow!(
            "acf-max-lag ({}) must be less than the number of m6A observations ({}).",
            max_lag,
            x.len()
        ));
    }

    let m = max_lag + 1;

    let len_x_usize = x.len();
    let len_x: T = From::from(len_x_usize as u32);
    let sum: T = From::from(0.0);

    let sum_x: T = x.iter().fold(sum, |sum, &xi| sum + xi);
    let mean_x: T = sum_x / len_x;

    let mut y: Vec<T> = vec![From::from(0.0); m];

    for t in 0..m {
        y[t] = x
            .into_par_iter()
            .enumerate()
            .take(len_x_usize - t)
            .map(|(i, xi)| {
                let xi = *xi - mean_x;
                let xi_t = x[i + t] - mean_x;
                (xi * xi_t) / len_x
            })
            .sum();
        // we need y[0] to calculate the correlations, so we set it to 1.0 at the end
        if !covariance && t > 0 {
            y[t] = y[t] / y[0];
        }
    }
    if !covariance {
        y[0] = From::from(1.0);
    }
    Ok(y)
}

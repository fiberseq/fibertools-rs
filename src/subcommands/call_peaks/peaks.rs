use super::fdr::FdrEntry;
use super::chrom_names_and_lengths;
use crate::cli::CallPeaksOptions;
use crate::subcommands::pileup::{FiberseqPileup, FiberseqPileupOptions, FireTrackOptions};
use crate::utils::bio_io;
use anyhow::Result;
use std::io::Write;

/// A peak representing a local maximum in FIRE scores
#[derive(Debug)]
pub struct Peak<'a> {
    pub chrom: String,
    /// Start position of the local maximum region (0-based, inclusive)
    /// This is the median start position of underlying FIRE elements
    pub start: usize,
    /// End position of the local maximum region (0-based, exclusive)
    /// This is the median end position of underlying FIRE elements
    pub end: usize,
    /// FIRE score at this position
    pub score: f32,
    /// FDR value at this position
    pub fdr: f64,
    /// Index in the pileup track where the peak was called (for retrieving FIRE elements)
    pub peak_index: usize,
    /// Reference to the pileup track containing FIRE elements
    pub pileup: &'a FiberseqPileup<'a>,
}

impl<'a> std::fmt::Display for Peak<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // BED6 format: chrom start end name score strand
        write!(
            f,
            "{}\t{}\t{}\t{:.6}\t{:.0}\t.",
            self.chrom, self.start, self.end, self.fdr, self.score
        )
    }
}

impl<'a> Peak<'a> {
    /// Find local maxima from a pileup
    /// For consecutive positions with identical rolling max scores, save the entire region
    /// Only keeps peaks with FDR <= 0.05
    pub fn from_pileup(pileup: &'a FiberseqPileup<'a>, fdr_table: &[FdrEntry]) -> Vec<Self> {
        let mut peaks = Vec::new();

        let scores = &pileup.all_data.scores;
        let rolling_max_scores = pileup.rolling_max.as_ref()
            .expect("Rolling max scores should be calculated");

        let mut consecutive_maxima = Vec::new();

        for i in 0..scores.len() {
            // Skip positions with no score or negative scores
            if scores[i] < 0.0 {
                // If we have consecutive maxima, create a peak for the region
                if !consecutive_maxima.is_empty() {
                    if let Some(peak) = Self::create_peak_if_significant(pileup, consecutive_maxima.as_slice(), fdr_table) {
                        peaks.push(peak);
                    }
                    consecutive_maxima.clear();
                }
                continue;
            }

            // Check if this position is a local maximum
            // (score equals the rolling max at this position)
            if (scores[i] - rolling_max_scores[i]).abs() < 1e-6 {
                consecutive_maxima.push(i);
            } else {
                // If we have consecutive maxima, create a peak for the region
                if !consecutive_maxima.is_empty() {
                    if let Some(peak) = Self::create_peak_if_significant(pileup, consecutive_maxima.as_slice(), fdr_table) {
                        peaks.push(peak);
                    }
                    consecutive_maxima.clear();
                }
            }
        }

        // Handle any remaining consecutive maxima at the end
        if !consecutive_maxima.is_empty() {
            if let Some(peak) = Self::create_peak_if_significant(pileup, consecutive_maxima.as_slice(), fdr_table) {
                peaks.push(peak);
            }
        }

        peaks
    }

    /// Create a peak from a region of consecutive maxima if FDR <= 0.05
    /// Returns None if the peak doesn't meet the FDR threshold
    /// Peak boundaries are determined by the median start/end of underlying FIRE elements
    fn create_peak_if_significant(pileup: &'a FiberseqPileup<'a>, positions: &[usize], fdr_table: &[FdrEntry]) -> Option<Self> {
        if positions.is_empty() {
            return None;
        }

        // Use the middle position to get the score and FDR
        let middle_idx = positions.len() / 2;
        let middle_pos = positions[middle_idx];
        let score = pileup.all_data.scores[middle_pos];
        let fdr = Self::lookup_fdr(score, fdr_table);

        // Only keep peaks with FDR <= 0.05
        if fdr <= 0.05 {
            // Collect all FIRE elements from the local max region
            let (start, end) = if let Some(ref fire_elements_vec) = pileup.all_data.fire_elements {
                let mut starts = Vec::new();
                let mut ends = Vec::new();

                // Collect FIRE elements from all positions in the local max region
                for &pos in positions {
                    for fire_elem in &fire_elements_vec[pos] {
                        starts.push(fire_elem.start);
                        ends.push(fire_elem.end);
                    }
                }

                if starts.is_empty() {
                    // Fallback: no FIRE elements found, use the local max region boundaries
                    let start = pileup.chrom_start + positions[0];
                    let end = pileup.chrom_start + positions[positions.len() - 1] + 1;
                    (start, end)
                } else {
                    // Calculate median start and end from FIRE elements
                    starts.sort_unstable();
                    ends.sort_unstable();
                    let median_start = starts[starts.len() / 2] as usize;
                    let median_end = ends[ends.len() / 2] as usize;
                    (median_start, median_end)
                }
            } else {
                // Fallback: FIRE element tracking not enabled
                let start = pileup.chrom_start + positions[0];
                let end = pileup.chrom_start + positions[positions.len() - 1] + 1;
                (start, end)
            };

            Some(Self {
                chrom: pileup.chrom.clone(),
                start,
                end,
                score,
                fdr,
                peak_index: middle_pos,
                pileup,
            })
        } else {
            None
        }
    }

    /// Look up FDR value for a given score
    fn lookup_fdr(score: f32, fdr_table: &[FdrEntry]) -> f64 {
        // Binary search to find the nearest threshold
        // We want the largest threshold that is <= score
        let idx = fdr_table.binary_search_by(|entry| {
            entry.threshold.partial_cmp(&(score as f64)).unwrap_or(std::cmp::Ordering::Equal)
        });

        match idx {
            Ok(i) => fdr_table[i].fdr,
            Err(i) => {
                if i == 0 {
                    fdr_table[0].fdr
                } else {
                    fdr_table[i - 1].fdr
                }
            }
        }
    }
}

/// Call peaks using the FDR table
///
/// This will:
/// 1. Run pileup on real data with FDR annotation
/// 2. Find local maxima
/// 3. Call peaks above FDR threshold
/// 4. Merge overlapping peaks
/// 5. Output results to BED file
pub fn call_peaks_with_fdr(
    opts: &mut CallPeaksOptions,
    bam: &mut rust_htslib::bam::IndexedReader,
    header: &rust_htslib::bam::HeaderView,
    fdr_table: &[FdrEntry],
) -> Result<()> {
    log::info!(
        "Calling peaks using FDR table with {} entries",
        fdr_table.len()
    );

    // Check if FDR table is empty
    if fdr_table.is_empty() {
        anyhow::bail!(
            "FDR table is empty. Cannot call peaks without an FDR table. \
            Please generate an FDR table first or provide a valid FDR table file."
        );
    }

    // Open output file and write header
    let mut writer = bio_io::writer(&opts.out)?;
    writeln!(writer, "#chrom\tstart\tend\tfdr\tscore\tstrand")?;

    let mut total_peaks = 0;

    for (chrom, chrom_len) in chrom_names_and_lengths(header)? {
        log::info!(
            "Finding peaks on chromosome {} with length {}",
            chrom,
            chrom_len
        );

        // Create FiberseqPileupOptions for peak calling
        let pileup_opts = FiberseqPileupOptions {
            fire_track_opts: FireTrackOptions {
                no_nuc: false,
                no_msp: false,
                m6a: false,
                cpg: false,
                fiber_coverage: true,
                shuffle: false,
                random_shuffle: false,
                shuffle_seed: None,
                rolling_max: Some(opts.window_size),
                track_fire_elements: true, // Enable FIRE element tracking for peak calling
            },
            rolling_max: Some(opts.window_size),
            haps: false,
            per_base: false,
            keep_zeros: false,
        };

        // Process fibers to build the track
        let mut pileup = FiberseqPileup::new(&chrom, 0, chrom_len as usize, pileup_opts, &None);
        let fibers = opts.input.fetch_fibers(bam, &chrom, None, None)?;
        pileup.add_fibers(fibers);

        // Find local maxima and filter by FDR <= 0.05
        let peaks = Peak::from_pileup(&pileup, fdr_table);
        log::info!(
            "Found {} significant peaks (FDR <= 0.05) on chromosome {}",
            peaks.len(),
            chrom
        );

        // Write peaks for this chromosome immediately
        for peak in &peaks {
            writeln!(writer, "{}", peak)?;
        }

        total_peaks += peaks.len();
    }

    log::info!(
        "Total peaks found across all chromosomes: {}",
        total_peaks
    );
    log::info!("Peaks written to {}", opts.out);

    Ok(())
}

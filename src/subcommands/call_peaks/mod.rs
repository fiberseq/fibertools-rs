mod fdr;
mod peaks;

pub use fdr::{
    fdr_table, read_fdr_table, write_fdr_table, FdrEntry, IncrementalFdrBuilder, PileupRecord,
};
pub use peaks::call_peaks;

use crate::cli::CallPeaksOptions;
use crate::fiber::FiberseqData;
use crate::subcommands::pileup::{FireTrack, FireTrackOptions};
use anyhow::{Context, Result};

pub fn run_call_peaks(opts: &mut CallPeaksOptions) -> Result<()> {
    log::info!("Starting FIRE peak calling");
    log::info!("  Input BAM: {}", opts.input.bam);
    log::info!("  Output: {}", opts.out);

    if let Some(min_frac) = opts.min_fire_frac {
        log::info!("  Using FIRE fraction mode: min_fire_frac = {}", min_frac);
    } else {
        log::info!("  Max FDR: {}", opts.max_fdr);
    }
    log::info!("  Window size: {}", opts.window_size);

    let mut bam = opts.input.indexed_bam_reader();
    let header = opts.input.header_view();

    // Generate or load FDR table (skip if using FIRE fraction mode)
    let fdr_table = if opts.min_fire_frac.is_some() {
        // FIRE fraction mode: use empty FDR table (won't be used for filtering)
        log::info!("  Skipping FDR calculation (using FIRE fraction threshold)");
        Vec::new()
    } else {
        // FDR mode: generate or load FDR table
        if let Some(ref shuffled) = opts.shuffled {
            log::info!("  Shuffled positions file: {}", shuffled);
        } else {
            log::info!("  Using random shuffling for FDR calculation");
        }

        if let Some(ref fdr_table_path) = opts.fdr_table {
            log::info!("Loading FDR table from: {}", fdr_table_path);
            read_fdr_table(fdr_table_path)?
        } else {
            // Generate pileup incrementally to avoid OOM
            log::info!("Running pileup for real and shuffled data (incremental mode)...");
            fdr_table(opts, &mut bam, &header)?
        }
    };

    // Write FDR table if requested (and if we generated one)
    if let Some(ref fdr_out) = opts.fdr_table_out {
        if !fdr_table.is_empty() {
            log::info!("Writing FDR table to: {}", fdr_out);
            write_fdr_table(&fdr_table, fdr_out)?;
        }
    }

    // Call peaks using the FDR table (or FIRE fraction)
    call_peaks(opts, &mut bam, &header, &fdr_table)?;

    log::info!("FIRE peak calling completed");
    Ok(())
}

/// Default minimum coverage threshold (matching Python's MIN_COVERAGE default)
const DEFAULT_MIN_COVERAGE: i32 = 4;

/// Get chromosome names and lengths from BAM header
///
/// # Returns
/// Vector of tuples (chromosome_name, chromosome_length)
pub fn chrom_names_and_lengths(
    header: &rust_htslib::bam::HeaderView,
) -> Result<Vec<(String, i64)>> {
    let mut chroms = Vec::new();
    for chrom in header.target_names() {
        let chrom_str = String::from_utf8_lossy(chrom).to_string();
        let tid = header.tid(chrom).context("Failed to get target ID")?;
        let chrom_len = header
            .target_len(tid)
            .context("Failed to get target length")? as i64;
        chroms.push((chrom_str, chrom_len));
    }
    Ok(chroms)
}

/// Process a single chromosome and return PileupRecords for both real and shuffled
/// Returns (real_records, shuffled_records) processed at the same positions
///
/// # Arguments
/// * `chrom` - Chromosome name
/// * `chrom_len` - Chromosome length
/// * `all_fibers` - Pre-loaded fibers from the chromosome
/// * `opts` - Call peaks options
fn process_chromosome_pileup_both(
    chrom: &str,
    chrom_len: i64,
    all_fibers: &[FiberseqData],
    opts: &CallPeaksOptions,
) -> Result<(Vec<PileupRecord>, Vec<PileupRecord>)> {
    if all_fibers.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }

    log::info!("Processing chromosome {} (length: {})", chrom, chrom_len);

    // Create fire track options - to calculate FIRE scores
    // Note: Python uses --no-msp --no-nuc flags in shuffled_pileup_chromosome rule
    // This means FIRE scores are calculated using only fiber_coverage, not MSP/NUC
    let real_opts = FireTrackOptions {
        no_nuc: true, // Match Python: --no-nuc
        no_msp: true, // Match Python: --no-msp
        m6a: false,
        cpg: false,
        fiber_coverage: true,
        shuffle: false,
        random_shuffle: false,
        shuffle_seed: None,
        rolling_max: None,
        track_fire_elements: false, // No tracking needed for FDR calculation
    };
    let mut real_track = FireTrack::new(chrom.to_string(), 0, chrom_len as usize, real_opts, &None);

    // Process all fibers and update real track
    // Note: FireTrack will automatically store fiber info in fibers_seen
    for fiber in all_fibers {
        real_track.update_with_fiber(fiber);
    }

    // Calculate scores for real fibers
    real_track.calculate_scores(Some(-1));

    // Calculate coverage thresholds based on median and standard deviation
    let (median, std_dev, pos_cov) = real_track.median_and_std_coverage();

    // Apply sd_cov thresholds if max_cov/min_cov are not explicitly set
    // Match Python behavior: minimum coverage defaults to 4
    let min_cov_threshold = opts.min_cov.unwrap_or_else(|| {
        let calculated_min = (median - opts.sd_cov * std_dev).round() as i32;
        calculated_min.max(DEFAULT_MIN_COVERAGE)
    });
    let max_cov_threshold = opts
        .max_cov
        .unwrap_or_else(|| (median + opts.sd_cov * std_dev).round() as i32);

    log::info!(
        "  Coverage: median={:.1}, std_dev={:.1} ({:.1} SDs), range=[{}, {}]",
        median,
        std_dev,
        opts.sd_cov,
        min_cov_threshold,
        max_cov_threshold
    );
    log::debug!("  Real: pos_with_cov={}", pos_cov);

    // Generate shuffled positions if not using a file-based shuffle
    // Pass coverage thresholds to avoid placing shuffled fibers in extreme coverage regions
    let generated_shuffle = Some(real_track.generate_shuffled_positions(
        Some(42),
        Some(min_cov_threshold),
        Some(max_cov_threshold),
    ));

    // Create shuffled track with shuffle enabled
    // Note: Must match real_opts to ensure consistent score calculation
    let shuffled_opts = FireTrackOptions {
        no_nuc: true, // Match Python: --no-nuc
        no_msp: true, // Match Python: --no-msp
        m6a: false,
        cpg: false,
        fiber_coverage: true,
        shuffle: true,
        random_shuffle: false, // We have explicit shuffle positions now
        shuffle_seed: None,
        rolling_max: None,
        track_fire_elements: false, // No tracking needed for FDR calculation
    };
    let mut shuffled_track = FireTrack::new(
        chrom.to_string(),
        0,
        chrom_len as usize,
        shuffled_opts,
        &generated_shuffle,
    );

    // Process the same fibers for shuffled track
    for fiber in all_fibers {
        shuffled_track.update_with_fiber(fiber);
    }

    // Calculate scores for shuffled track
    shuffled_track.calculate_scores(Some(-1));

    // Log median coverage statistics for shuffled
    let (shuffled_median, shuffled_pos_cov) = shuffled_track.median_coverage();
    log::info!(
        "  Shuffled: median_cov={:.1}, pos_with_cov={}",
        shuffled_median,
        shuffled_pos_cov
    );

    // Extract pileup records at positions where EITHER real or shuffled has coverage
    // This ensures both tracks use the same position set
    let mut real_records = Vec::new();
    let mut shuffled_records = Vec::new();

    for i in 0..real_track.track_len {
        let real_cov = real_track.coverage[i];
        let shuffled_cov = shuffled_track.coverage[i];
        let real_score = real_track.scores[i];
        let shuffled_score = shuffled_track.scores[i];

        // Skip positions with:
        // - No coverage (cov <= 0)
        // - Invalid scores (score < 0), which indicates low FIRE coverage
        // - Coverage outside min/max thresholds (matching Python filtering)
        let skip_real = real_cov <= 0
            || real_score < 0.0
            || real_cov < min_cov_threshold
            || real_cov > max_cov_threshold;
        let skip_shuffled = shuffled_cov <= 0
            || shuffled_score < 0.0
            || shuffled_cov < min_cov_threshold
            || shuffled_cov > max_cov_threshold;

        // Add real record at positions with valid coverage and score
        if !skip_real {
            real_records.push(PileupRecord {
                start: i as u64,
                end: (i + 1) as u64,
                coverage: real_cov as u32,
                fire_coverage: real_track.fire_coverage[i] as u32,
                score: real_score as f64,
            });
        }

        // Add shuffled record at positions with valid coverage and score
        if !skip_shuffled {
            shuffled_records.push(PileupRecord {
                start: i as u64,
                end: (i + 1) as u64,
                coverage: shuffled_cov as u32,
                fire_coverage: shuffled_track.fire_coverage[i] as u32,
                score: shuffled_score as f64,
            });
        }
    }

    // log the lowest and highest scores in real and shuffled tracks
    for (label, track) in &[("real", &real_records), ("shuffle", &shuffled_records)] {
        let n_neg_one = track.iter().filter(|x| x.score < 0.0).count();
        let min = track
            .iter()
            .filter(|x| x.score > 0.0)
            .fold(f64::INFINITY, |a, b| a.min(b.score));
        let max = track.iter().fold(f64::NEG_INFINITY, |a, b| a.max(b.score));
        log::debug!(
            "    Track {}: score range = [{:.3}, {:.3}], n_scores<0={}",
            label,
            min,
            max,
            n_neg_one
        );
    }

    Ok((real_records, shuffled_records))
}

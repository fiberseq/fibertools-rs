use crate::cli::CallPeaksOptions;
use anyhow::{Context, Result};
use itertools::Itertools;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::io::Write;

/// FDR table entry mapping FIRE scores to FDR values
#[derive(Debug, Clone)]
pub struct FdrEntry {
    pub threshold: f64,
    pub fdr: f64,
    pub shuffled_bp: f64,
    pub real_bp: f64,
}

/// Pileup record structure (simplified for now)
#[derive(Debug, Clone)]
pub struct PileupRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub coverage: u32,
    pub fire_coverage: u32,
    pub score: f64,
}

/// Calculate FDR from aggregated FIRE scores
/// This follows the Python logic in fdr_from_fire_scores()
fn fdr_from_fire_scores(
    fire_scores: &[(f64, bool, u64)],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut vs = Vec::new(); // shuffled bp
    let mut rs = Vec::new(); // real bp
    let mut ts = Vec::new(); // thresholds
    let mut cur_r = 0.0;
    let mut cur_v = 0.0;
    let mut pre_score = -1.0;
    let mut first = true;

    for &(score, is_real, bp) in fire_scores {
        // save the counts and thresholds as long as we have counts
        if score != pre_score && cur_r > 0.0 && !first {
            rs.push(cur_r);
            vs.push(cur_v);
            ts.push(pre_score);
        }

        // don't add negative scores to the fdr data, since they have no coverage
        if score < 0.0 {
            break;
        }

        // update the counts
        if is_real {
            cur_r += bp as f64;
        } else {
            cur_v += bp as f64;
        }

        // prepare for next iteration
        pre_score = score;
        first = false;
    }

    // add the last threshold with an FDR of 1
    rs.push(1.0);
    vs.push(1.0);
    ts.push(-1.0);

    // calculate FDRs
    let fdrs: Vec<f64> = vs
        .iter()
        .zip(rs.iter())
        .map(|(&v, &r)| {
            let fdr = v / r;
            if fdr > 1.0 {
                1.0
            } else {
                fdr
            }
        })
        .collect();

    (ts, fdrs, vs, rs)
}

/// Create FDR table from fire scores
/// This follows the Python logic in fdr_table_from_scores()
fn fdr_table_from_scores(fire_scores: &[(f64, bool, u64)]) -> Vec<FdrEntry> {
    let (thresholds, fdrs, shuffled_bps, real_bps) = fdr_from_fire_scores(fire_scores);

    let mut entries: Vec<FdrEntry> = thresholds
        .iter()
        .zip(fdrs.iter())
        .zip(shuffled_bps.iter())
        .zip(real_bps.iter())
        .map(|(((&threshold, &fdr), &shuffled_bp), &real_bp)| FdrEntry {
            threshold,
            fdr,
            shuffled_bp,
            real_bp,
        })
        .collect();

    // simplify the results - group by FDR and keep last
    entries.sort_by(|a, b| a.fdr.partial_cmp(&b.fdr).unwrap());
    deduplicate_by_key(&mut entries, |e| (e.fdr * 1000000.0) as i64);

    // group by shuffled_bp and keep last
    entries.sort_by(|a, b| a.shuffled_bp.partial_cmp(&b.shuffled_bp).unwrap());
    deduplicate_by_key(&mut entries, |e| e.shuffled_bp as i64);

    // group by real_bp and keep last
    entries.sort_by(|a, b| a.real_bp.partial_cmp(&b.real_bp).unwrap());
    deduplicate_by_key(&mut entries, |e| e.real_bp as i64);

    // round thresholds to 2 decimal places
    for entry in &mut entries {
        entry.threshold = (entry.threshold * 100.0).round() / 100.0;
    }

    // group by threshold and keep last
    entries.sort_by(|a, b| a.threshold.partial_cmp(&b.threshold).unwrap());
    deduplicate_by_key(&mut entries, |e| (e.threshold * 100.0) as i64);

    // sort by threshold ascending (needed for binary search later)
    entries.sort_by(|a, b| a.threshold.partial_cmp(&b.threshold).unwrap());

    log::info!("FDR table has {} entries", entries.len());
    if !entries.is_empty() {
        log::debug!(
            "First FDR entry: threshold={:.2}, FDR={:.4}",
            entries[0].threshold,
            entries[0].fdr
        );
        log::debug!(
            "Last FDR entry: threshold={:.2}, FDR={:.4}",
            entries.last().unwrap().threshold,
            entries.last().unwrap().fdr
        );
    }

    entries
}

/// Helper function to deduplicate entries by key, keeping the last occurrence
fn deduplicate_by_key<T, K: Eq + std::hash::Hash>(entries: &mut Vec<T>, key_fn: impl Fn(&T) -> K)
where
    T: Clone,
{
    let mut seen = HashMap::new();
    for (idx, entry) in entries.iter().enumerate() {
        seen.insert(key_fn(entry), idx);
    }
    let mut keep_indices: Vec<_> = seen.values().copied().collect();
    keep_indices.sort_unstable();

    let mut result = Vec::with_capacity(keep_indices.len());
    for &idx in &keep_indices {
        result.push(entries[idx].clone());
    }
    *entries = result;
}

/// Aggregate pileup data by FIRE score
/// Returns HashMap mapping score to total base pairs
/// Uses full float precision for grouping (matching Python behavior)
fn aggregate_pileup_by_score(
    pileup_records: &[PileupRecord],
    max_cov: Option<u32>,
    min_cov: Option<u32>,
) -> HashMap<ordered_float::NotNan<f64>, u64> {
    use ordered_float::NotNan;

    let mut score_counts: HashMap<NotNan<f64>, u64> = HashMap::new();

    for record in pileup_records {
        // filter by coverage if specified
        if let Some(max) = max_cov {
            if record.coverage > max {
                continue;
            }
        }
        if let Some(min) = min_cov {
            if record.coverage < min {
                continue;
            }
        }

        let bp = record.end - record.start;
        // Use full precision, don't round yet (matches Python line 178)
        if let Ok(score_key) = NotNan::new(record.score) {
            *score_counts.entry(score_key).or_insert(0) += bp;
        }
    }

    score_counts
}

/// Make FDR table from real and shuffled pileup data
fn make_fdr_table(opts: &mut CallPeaksOptions) -> Result<Vec<FdrEntry>> {
    log::info!("Generating FDR table from pileup data");

    // First pass: generate pileup for real data
    log::info!("Running pileup on real data...");
    let real_pileup = run_pileup_for_peaks(opts, false)?;
    log::info!("Real pileup generated {} records", real_pileup.len());

    // Second pass: generate pileup for shuffled data
    // Will use file-based shuffle if provided, otherwise random shuffle
    log::info!("Running pileup on shuffled data...");
    let shuffled_pileup = run_pileup_for_peaks(opts, true)?;
    log::info!(
        "Shuffled pileup generated {} records",
        shuffled_pileup.len()
    );

    // Aggregate by score
    let real_scores = aggregate_pileup_by_score(&real_pileup, opts.max_cov, opts.min_cov);
    let shuffled_scores = aggregate_pileup_by_score(&shuffled_pileup, opts.max_cov, opts.min_cov);

    // Combine and sort by score descending
    let mut fire_scores: Vec<(f64, bool, u64)> = Vec::new();
    for (score_notnan, bp) in real_scores {
        fire_scores.push((score_notnan.into_inner(), true, bp));
    }
    for (score_notnan, bp) in shuffled_scores {
        fire_scores.push((score_notnan.into_inner(), false, bp));
    }
    fire_scores.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap()); // descending by score

    // Calculate sums for logging
    let real_mbp: f64 = fire_scores
        .iter()
        .filter(|(_, is_real, _)| *is_real)
        .map(|(_, _, bp)| *bp as f64)
        .sum::<f64>()
        / 1_000_000.0;
    let shuffled_mbp: f64 = fire_scores
        .iter()
        .filter(|(_, is_real, _)| !*is_real)
        .map(|(_, _, bp)| *bp as f64)
        .sum::<f64>()
        / 1_000_000.0;
    log::info!("Real data: {:.2} Mbp", real_mbp);
    log::info!("Shuffled data: {:.2} Mbp", shuffled_mbp);

    // Create FDR table
    let fdr_table = fdr_table_from_scores(&fire_scores);

    // Check if we have any thresholds below max_fdr
    if let Some(min_fdr_entry) = fdr_table
        .iter()
        .min_by(|a, b| a.fdr.partial_cmp(&b.fdr).unwrap())
    {
        if min_fdr_entry.fdr > opts.max_fdr {
            anyhow::bail!(
                "No FIRE score threshold has an FDR < {}. Check the input Fiber-seq data with the QC pipeline and make sure you are using WGS Fiber-seq data.",
                opts.max_fdr
            );
        }
    }

    Ok(fdr_table)
}

/// Process a single chromosome and return PileupRecords
/// Processes the entire chromosome at once (no windowing)
fn process_chromosome_pileup(
    chrom: &str,
    chrom_len: i64,
    bam: &mut rust_htslib::bam::IndexedReader,
    opts: &CallPeaksOptions,
    shuffled_fibers: &Option<crate::subcommands::pileup::ShuffledFibers>,
    shuffled: bool,
) -> Result<Vec<PileupRecord>> {
    use crate::fiber::FiberseqData;
    use crate::subcommands::pileup::{FireTrack, FireTrackOptions};

    log::debug!("Processing chromosome {} (length: {})", chrom, chrom_len);

    // Check if region has data
    bam.fetch((chrom, 0, chrom_len))?;
    let mut tmp_records = bam.records();
    if tmp_records.next().is_none() {
        return Ok(Vec::new());
    }

    // Fetch the data again for processing
    bam.fetch((chrom, 0, chrom_len))?;
    let records = bam.records();

    // Create fire track for entire chromosome
    let shuffled_ref = if shuffled { shuffled_fibers } else { &None };

    let fire_track_opts = FireTrackOptions {
        no_nuc: !opts.include_nuc_msp,
        no_msp: !opts.include_nuc_msp,
        m6a: false,
        cpg: false,
        fiber_coverage: false,
        shuffle: shuffled,
        random_shuffle: shuffled && shuffled_fibers.is_none(), // Use random shuffle if no file provided
        shuffle_seed: Some(42),                                // TODO: Could make this a CLI option
        rolling_max: None,
    };

    let mut fire_track = FireTrack::new(0, chrom_len as usize, fire_track_opts, shuffled_ref);

    // Process records
    opts.input
        .filters
        .filter_on_bit_flags(records)
        .chunks(1000)
        .into_iter()
        .for_each(|chunk| {
            let chunk: Vec<_> = chunk.collect();
            let fibers: Vec<FiberseqData> =
                FiberseqData::from_records(chunk, &opts.input.header_view(), &opts.input.filters);

            for fiber in fibers {
                // Skip if shuffled and fiber not in shuffled set
                if shuffled {
                    if let Some(ref sf) = shuffled_fibers {
                        if !sf.has_fiber(&fiber) {
                            return;
                        }
                    }
                }

                fire_track.update_with_fiber(&fiber);
            }
        });

    // Calculate scores
    fire_track.calculate_scores();

    // Extract pileup records
    let mut records = Vec::new();
    for i in 0..fire_track.track_len {
        let score = fire_track.scores[i];
        let coverage = fire_track.coverage[i];
        let fire_coverage = fire_track.fire_coverage[i];

        // Apply coverage filters if specified
        if let Some(max_cov) = opts.max_cov {
            if coverage as u32 > max_cov {
                continue;
            }
        }
        if let Some(min_cov) = opts.min_cov {
            if (coverage as u32) < min_cov {
                continue;
            }
        }

        // Only include positions with valid scores or coverage
        if score >= 0.0 || coverage > 0 {
            records.push(PileupRecord {
                chrom: chrom.to_string(),
                start: i as u64,
                end: (i + 1) as u64,
                coverage: coverage as u32,
                fire_coverage: fire_coverage as u32,
                score: score as f64,
            });
        }
    }

    Ok(records)
}

/// Run pileup and yield records per chromosome
/// Processes each chromosome completely before moving to the next
fn run_pileup_for_peaks(opts: &mut CallPeaksOptions, shuffled: bool) -> Result<Vec<PileupRecord>> {
    use crate::subcommands::pileup::ShuffledFibers;

    log::info!(
        "Running pileup for {} data",
        if shuffled { "shuffled" } else { "real" }
    );

    let mut bam = opts.input.indexed_bam_reader();
    let header = opts.input.header_view();

    // Load shuffled fibers if file is provided
    // If shuffled=true but no file, will use random shuffling instead
    let shuffled_fibers = match &opts.shuffled {
        Some(path) => {
            log::info!("Loading shuffled positions from file: {}", path);
            Some(ShuffledFibers::new(path)?)
        }
        None if shuffled => {
            log::info!("Using random shuffling (no shuffle file provided)");
            None // Will trigger random_shuffle in FireTrackOptions
        }
        None => None,
    };

    let mut all_records = Vec::new();

    // Process each chromosome
    for chrom in header.target_names() {
        let chrom_str = String::from_utf8_lossy(chrom).to_string();
        let tid = bam.header().tid(chrom).context("Failed to get target ID")?;
        let chrom_len = bam
            .header()
            .target_len(tid)
            .context("Failed to get target length")? as i64;

        let chrom_records = process_chromosome_pileup(
            &chrom_str,
            chrom_len,
            &mut bam,
            opts,
            &shuffled_fibers,
            shuffled,
        )?;

        log::debug!(
            "Chromosome {} yielded {} records",
            chrom_str,
            chrom_records.len()
        );
        all_records.extend(chrom_records);
    }

    log::info!("Collected {} pileup records total", all_records.len());
    Ok(all_records)
}

pub fn run_call_peaks(opts: &mut CallPeaksOptions) -> Result<()> {
    log::info!("Starting FIRE peak calling");
    log::info!("  Input BAM: {}", opts.input.bam);
    log::info!("  Output: {}", opts.out);
    log::info!("  Max FDR: {}", opts.max_fdr);
    log::info!("  Window size: {}", opts.window_size);

    if let Some(ref shuffled) = opts.shuffled {
        log::info!("  Shuffled positions file: {}", shuffled);
    } else {
        log::info!("  Using random shuffling for FDR calculation");
    }

    // Generate or load FDR table
    let fdr_table = if let Some(ref fdr_table_path) = opts.fdr_table {
        log::info!("Loading FDR table from: {}", fdr_table_path);
        // TODO: load from file
        Vec::new()
    } else {
        // Always generate FDR table (will use random shuffle if no file provided)
        make_fdr_table(opts)?
    };

    // Write FDR table if requested
    if let Some(ref fdr_out) = opts.fdr_table_out {
        log::info!("Writing FDR table to: {}", fdr_out);
        write_fdr_table(&fdr_table, fdr_out)?;
    }

    // TODO: Implement peak calling logic
    // 1. Run pileup on real data with FDR annotation
    // 2. Find local maxima
    // 3. Call peaks
    // 4. Merge peaks
    // 5. Output results

    log::info!("FIRE peak calling completed");
    Ok(())
}

/// Write FDR table to TSV file
fn write_fdr_table(fdr_table: &[FdrEntry], path: &str) -> Result<()> {
    let mut writer =
        crate::utils::bio_io::writer(path).context("Failed to create FDR table output file")?;

    writeln!(writer, "threshold\tFDR\tshuffled_bp\treal_bp")?;
    for entry in fdr_table {
        writeln!(
            writer,
            "{:.2}\t{:.6}\t{:.0}\t{:.0}",
            entry.threshold, entry.fdr, entry.shuffled_bp, entry.real_bp
        )?;
    }

    Ok(())
}

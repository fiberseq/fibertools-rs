/// This module is used to extract the fire calls as well as nucs and msps from a bam file
/// for every position in the bam file and output the results to a bed file.
/// all calculations are done in total as well as for haplotype 1 and haplotype 2.
use self::fiber::FiberseqData;
use self::fiber::FiberseqRecords;
use super::cli::FireTrackOptions;
use super::*;
use anyhow;
//use polars::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;

const MIN_FIRE_COVERAGE: i32 = 4;
const MIN_FIRE_QUAL: u8 = 229; // floor(255*0.9)

pub struct FireTrack {
    pub chrom_len: usize,
    pub scores: Vec<f32>,
    pub coverage: Vec<i32>,
    pub fire_coverage: Vec<i32>,
}

impl FireTrack {
    pub fn new(fibers: &[FiberseqData], chrom_len: usize, shuffle_offsets: &[i64]) -> Self {
        let scores = vec![-1.0; chrom_len];
        let coverage = vec![0; chrom_len];
        let fire_coverage = vec![0; chrom_len];
        let mut chrom_data = Self {
            chrom_len,
            scores,
            coverage,
            fire_coverage,
        };
        chrom_data.update_with_fibers(fibers, shuffle_offsets);
        /*
        let df = df!(
            "coverage" => &chrom_data.coverage,
            "fire_coverage" => &chrom_data.fire_coverage,
            "scores" => &chrom_data.scores
        );
         */
        chrom_data
    }

    fn update_with_fibers(&mut self, fibers: &[FiberseqData], shuffle_offsets: &[i64]) {
        for (fiber, shuffle_offset) in fibers.iter().zip(shuffle_offsets) {
            // calculate the coverage
            for i in fiber.record.reference_start()..fiber.record.reference_end() {
                let i = (i + shuffle_offset) as usize;
                self.coverage[i] += 1;
            }
            // calculate the fire coverage and fire score
            for (_, _, _, q, r) in &fiber.msp {
                match r {
                    Some((rs, re, _)) => {
                        if q < MIN_FIRE_QUAL {
                            continue;
                        }
                        let score_update = (1.0 / (q as f32)).log10();
                        for i in rs..re {
                            let i = (i + shuffle_offset) as usize;
                            self.fire_coverage[i] += 1;
                            self.scores[i] += score_update;
                        }
                    }
                    _ => continue,
                }
            }
        }

        // correct the score for coverage and normalize to 0-100
        for i in 0..self.chrom_len {
            if self.fire_coverage[i] < MIN_FIRE_COVERAGE {
                self.scores[i] = -1.0;
            } else {
                self.scores[i] = -50.0 * self.scores[i] / self.coverage[i] as f32;
            }
        }
    }

    pub fn header() -> String {
        format!(
            "{}{}{}{}{}{}",
            "#chrom\tstart\tend\t",
            "fire_coverage\tcoverage\t",
            "score\tFDR\tlog_FDR\t",
            "fire_coverage_H1\tcoverage_H1\tscore_H1\tFDR_H1\tlog_FDR_H1\t",
            "fire_coverage_H2\tcoverage_H2\tscore_H2\tFDR_H2\tlog_FDR_H2\t",
            "max_window_score\tis_local_max\n"
        )
    }
}

/// extract existing fire calls into a bed9+ like file
pub fn fire_track(fire_track_opts: &FireTrackOptions) -> Result<(), anyhow::Error> {
    let mut bam = bio_io::bam_reader(&fire_track_opts.bam, 8);
    let mut out_buffer = bio_io::writer(&fire_track_opts.out)?;
    // add the header
    out_buffer.write_all(FireTrack::header().as_bytes())?;

    for rec in FiberseqRecords::new(&mut bam, 0) {
        if rec.record.is_secondary() || rec.record.is_supplementary() || rec.record.is_unmapped() {
            continue;
        }
    }
    todo!();
}

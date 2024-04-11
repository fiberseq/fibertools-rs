use self::cli::PileupOptions;
use self::fiber::FiberseqData;
use super::*;
use anyhow::{self, Ok};
use rayon::prelude::*;
//use polars::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{FetchDefinition, IndexedReader};

/// This module is used to extract the fire calls as well as nucs and msps from a bam file
/// for every position in the bam file and output the results to a bed file.
/// all calculations are done in total as well as for haplotype 1 and haplotype 2.

const MIN_FIRE_COVERAGE: i32 = 4;
const MIN_FIRE_QUAL: u8 = 229; // floor(255*0.9)

pub struct FireTrack {
    pub chrom_len: usize,
    pub raw_scores: Vec<f32>,
    pub scores: Vec<f32>,
    pub coverage: Vec<i32>,
    pub fire_coverage: Vec<i32>,
    pub msp_coverage: Vec<i32>,
    pub nuc_coverage: Vec<i32>,
    pub cpg_coverage: Vec<i32>,
    pub m6a_coverage: Vec<i32>,
}

impl FireTrack {
    pub fn new(chrom_len: usize) -> Self {
        let raw_scores = vec![-1.0; chrom_len];
        let scores = vec![-1.0; chrom_len];
        Self {
            chrom_len,
            raw_scores,
            scores,
            coverage: vec![0; chrom_len],
            fire_coverage: vec![0; chrom_len],
            msp_coverage: vec![0; chrom_len],
            nuc_coverage: vec![0; chrom_len],
            cpg_coverage: vec![0; chrom_len],
            m6a_coverage: vec![0; chrom_len],
        }
    }

    #[inline]
    fn add_range_set(array: &mut [i32], ranges: &super::bamranges::Ranges) {
        let shuffle_offset = 0;
        for (_, _, _, _, r) in ranges {
            match r {
                Some((rs, re, _)) => {
                    let re = if rs == re { re + 1 } else { re };
                    for i in rs..re {
                        let i = (i + shuffle_offset) as usize;
                        array[i] += 1;
                    }
                }
                _ => continue,
            }
        }
    }

    /// inline this function
    #[inline]
    fn update_with_fiber(&mut self, fiber: &FiberseqData) {
        let shuffle_offset = 0;
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
                    let score_update = (1.0 - q as f32 / 255.0).log10() * -50.0;
                    for i in rs..re {
                        let i = (i + shuffle_offset) as usize;
                        self.fire_coverage[i] += 1;
                        self.raw_scores[i] += score_update;
                    }
                }
                _ => continue,
            }
        }
        Self::add_range_set(&mut self.nuc_coverage, &fiber.nuc);
        Self::add_range_set(&mut self.msp_coverage, &fiber.msp);
        Self::add_range_set(&mut self.m6a_coverage, &fiber.m6a);
        Self::add_range_set(&mut self.cpg_coverage, &fiber.cpg);
    }

    pub fn calculate_scores(&mut self) {
        for i in 0..self.chrom_len {
            if self.fire_coverage[i] < MIN_FIRE_COVERAGE {
                self.scores[i] = -1.0;
            } else {
                self.scores[i] = self.raw_scores[i] / self.coverage[i] as f32;
            }
        }
    }
}

pub struct FiberseqPileup<'a> {
    pub all_data: FireTrack,
    pub hap1_data: FireTrack,
    pub hap2_data: FireTrack,
    pub chrom: String,
    pub chrom_len: usize,
    has_data: bool,
    pileup_opts: &'a PileupOptions,
}

impl<'a> FiberseqPileup<'a> {
    pub fn new(chrom: &str, chrom_len: usize, pileup_opts: &'a PileupOptions) -> Self {
        let all_data = FireTrack::new(chrom_len);
        let hap1_data = FireTrack::new(chrom_len);
        let hap2_data = FireTrack::new(chrom_len);
        Self {
            all_data,
            hap1_data,
            hap2_data,
            chrom: chrom.to_string(),
            chrom_len,
            has_data: false,
            pileup_opts,
        }
    }

    pub fn has_data(&self) -> bool {
        self.has_data
    }

    pub fn add_records(
        &mut self,
        records: bam::Records<'a, IndexedReader>,
    ) -> Result<(), anyhow::Error> {
        records
            .map(|r| r.unwrap())
            .filter(|r| !(r.is_secondary() || r.is_supplementary() || r.is_unmapped()))
            .chunks(1000)
            .into_iter()
            .map(|r| r.collect::<Vec<_>>())
            .for_each(|r| {
                let fibers: Vec<FiberseqData> = r
                    .par_iter()
                    .map(|r| FiberseqData::new(r.clone(), None, 0))
                    .collect();
                if !fibers.is_empty() {
                    self.has_data = true;
                }
                for fiber in fibers {
                    self.all_data.update_with_fiber(&fiber);
                    if fiber.get_hp() == "H1" {
                        self.hap1_data.update_with_fiber(&fiber);
                    } else if fiber.get_hp() == "H2" {
                        self.hap2_data.update_with_fiber(&fiber);
                    }
                }
            });
        Ok(())
    }

    pub fn header() -> String {
        format!(
            "{}\t{}\t{}\t{}\n",
            "#chrom\tpos",
            "coverage\tfire_coverage\tscore\tnuc_coverage\tmsp_coverage\tm6a_coverage\tcpg_coverage",
            "coverage_H1\tfire_coverage_H1\tscore_H1\tnuc_coverage_H1\tmsp_coverage_H1\tm6a_coverage_H1\tcpg_coverage_H1",
            "coverage_H2\tfire_coverage_H2\tscore_H2\tnuc_coverage_H2\tmsp_coverage_H2\tm6a_coverage_H2\tcpg_coverage_H2",
        )
    }

    pub fn calculate_scores(&mut self) {
        self.all_data.calculate_scores();
        self.hap1_data.calculate_scores();
        self.hap2_data.calculate_scores();
    }

    pub fn write(&self, out: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        if !self.has_data {
            return Ok(());
        }

        for i in 0..self.chrom_len {
            if !self.pileup_opts.keep_zeros && self.all_data.coverage[i] == 0 {
                continue;
            }
            let mut line = format!("{}\t{}\t", self.chrom, i);
            for data in [&self.all_data, &self.hap1_data, &self.hap2_data] {
                line += &format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
                    data.coverage[i],
                    data.fire_coverage[i],
                    data.scores[i],
                    data.nuc_coverage[i],
                    data.msp_coverage[i],
                    data.m6a_coverage[i],
                    data.cpg_coverage[i],
                );
            }
            // remove the last tab
            line.pop();
            line += "\n";

            out.write_all(line.as_bytes())?;
        }
        out.flush()?;
        Ok(())
    }
}

fn run_rgn(
    chrom: &str,
    rgn: FetchDefinition,
    bam: &mut IndexedReader,
    out: &mut Box<dyn Write>,
    pileup_opts: &PileupOptions,
) -> Result<(), anyhow::Error> {
    log::info!("Processing region: {}", chrom);
    let tid = bam.header().tid(chrom.as_bytes()).unwrap();
    let chrom_len = bam.header().target_len(tid).unwrap() as usize;
    bam.fetch(rgn)?;
    let records = bam.records();
    // make the pileup
    let mut pileup = FiberseqPileup::new(chrom, chrom_len, pileup_opts);
    pileup.add_records(records)?;
    pileup.calculate_scores();
    pileup.write(out)?;
    Ok(())
}

/// extract existing fire calls into a bed9+ like file
pub fn pileup_track(pileup_opts: &PileupOptions) -> Result<(), anyhow::Error> {
    // read in the bam from stdin or from a file
    let mut bam = bam::IndexedReader::from_path(pileup_opts.bam.clone())?;
    let header = bam.header().clone();
    bam.set_threads(8).unwrap();

    let mut out = bio_io::writer(&pileup_opts.out)?;
    // add the header
    out.write_all(FiberseqPileup::header().as_bytes())?;

    match &pileup_opts.rgn {
        // if a region is specified, only process that region
        Some(rgn) => {
            let (rgn, chrom) = region_parser(rgn);
            run_rgn(&chrom, rgn, &mut bam, &mut out, pileup_opts)?;
        }
        // if no region is specified, process all regions
        None => {
            for chrom in header.target_names() {
                let rgn = FetchDefinition::String(chrom);
                run_rgn(
                    &String::from_utf8_lossy(chrom),
                    rgn,
                    &mut bam,
                    &mut out,
                    pileup_opts,
                )?;
            }
        }
    }
    Ok(())
}

use crate::center::CenterPosition;
use crate::fiber;

use super::bam_writer;
use super::bio_io;
use super::cli::FootprintOptions;
use super::fiber::*;
use anyhow;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::*};
use serde::{Deserialize, Serialize};
use serde_yaml;
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct FootprintYaml {
    modules: Vec<(usize, usize)>,
}

impl FootprintYaml {
    pub fn check_for_valid_input(&self) -> Result<(), anyhow::Error> {
        if self.modules.is_empty() || self.modules.len() > 8 {
            return Err(anyhow::anyhow!(
                "YAML SPECIFICATION ERROR: Number of modules must be between 1 and 8: {self:?}"
            ));
        }
        let mut last_end = 0;
        for (start, end) in &self.modules {
            if *start < last_end {
                return Err(anyhow::anyhow!(
                    "YAML SPECIFICATION ERROR: modules are not sorted: {self:?}"
                ));
            }
            if *start != last_end {
                return Err(anyhow::anyhow!(
                    "YAML SPECIFICATION ERROR: modules are not contiguous or do not start at 0: {self:?}"
                ));
            }
            if *start > *end {
                return Err(anyhow::anyhow!(
                    "YAML SPECIFICATION ERROR: start > end: {self:?}"
                ));
            }
            last_end = *end;
        }
        Ok(())
    }

    pub fn max_pos(&self) -> usize {
        self.modules[self.modules.len() - 1].1
    }
}

pub struct ReferenceMotif<'a> {
    chrom: String,
    start: i64,
    end: i64,
    footprint: &'a FootprintYaml,
}

impl ReferenceMotif<'_> {
    pub fn spans(&self, start: i64, end: i64) -> bool {
        // check if coordinates span the motif
        self.start > start && self.end < end
    }
}

pub struct Footprint<'a> {
    n_spanning_fibers: usize,
    n_spanning_msps: usize,
    has_spanning_msp: Vec<bool>,
    footprint_codes: Vec<u16>,
    motif: &'a ReferenceMotif<'a>,
    fibers: Vec<&'a FiberseqData>,
}

impl<'a> Footprint<'a> {
    pub fn new(motif: &'a ReferenceMotif, in_fibers: &'a Vec<FiberseqData>) -> Self {
        let mut fibers = vec![];
        for fiber in in_fibers {
            // add if fiber spans the footprint
            if motif.spans(fiber.record.reference_start(), fiber.record.reference_end()) {
                fibers.push(fiber);
            }
        }

        let mut footprint = Self {
            n_spanning_fibers: fibers.len(),
            n_spanning_msps: 0,
            has_spanning_msp: vec![],
            footprint_codes: vec![],
            motif,
            fibers,
        };
        // add the number of msps that span the footprint
        footprint.spanning_msps();
        footprint.establish_footprints();
        footprint
    }

    fn spanning_msps(&mut self) {
        for fiber in self.fibers.iter() {
            let mut has_spanning_msp = false;
            for msp in &fiber.msp {
                // skip if there is no mapping of the msp
                match msp.3 {
                    Some((rs, re, _rl)) => {
                        if self.motif.spans(rs, re) {
                            self.n_spanning_msps += 1;
                            has_spanning_msp = true;
                            break;
                        }
                    }
                    None => continue,
                }
            }
            self.has_spanning_msp.push(has_spanning_msp);
        }
    }

    fn establish_footprints(&mut self) {
        for (fiber, has_msp) in self.fibers.iter().zip(self.has_spanning_msp.iter()) {
            let mut footprint_code = 0;
            if *has_msp {
                footprint_code = self.footprint_code(fiber);
            }
            self.footprint_codes.push(footprint_code);
        }
    }

    fn footprint_code(&self, fiber: &FiberseqData) -> u16 {
        // denote that we have a spanning msp
        let mut footprint_code = 1 << 0;
        let motif_end = self.motif.footprint.max_pos();

        // make a binary vector over the motif indicating the presence of an m6a
        let mut m6a_vec = vec![false; motif_end + 1];
        for m6a in fiber.m6a.reference_starts.iter().flatten() {
            if m6a < &self.motif.start || m6a > &self.motif.end {
                continue;
            }
            let m6a_rel_pos = m6a - self.motif.start;
            m6a_vec[m6a_rel_pos as usize] = true;
        }

        // mark modules within the motif footprint based on the presence of an m6a
        for (i, (start, end)) in self.motif.footprint.modules.iter().enumerate() {
            let count = m6a_vec[*start..*end].iter().filter(|&&x| x).count();
            // the module doesn't have an m6a, so it is a footprint
            if count == 0 {
                footprint_code |= 1 << (i + 1);
            }
        }

        footprint_code
    }
}

pub fn define_footprint(fiber: FiberseqData, bed_rec: CenterPosition, _modules: &FootprintYaml) {
    // check for an overlap
    if fiber.target_name != bed_rec.chrom
        || fiber.record.reference_start() > bed_rec.position
        || fiber.record.reference_end() < bed_rec.position
    {}
}

pub fn start_finding_footprints(opts: &FootprintOptions) -> Result<(), anyhow::Error> {
    let yaml_buff = bio_io::buffer_from(&opts.yaml)?;
    let yaml: FootprintYaml = serde_yaml::from_reader(yaml_buff)?;
    yaml.check_for_valid_input()?;
    log::debug!("YAML: {:?}", yaml);

    let bed_records = super::center::read_center_positions(&opts.bed)?;
    let bam = bio_io::bam_reader(&opts.bam, 1);
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let mut out = bam_writer(&opts.out, &bam, 8);
    let mut bam = rust_htslib::bam::IndexedReader::from_path(&opts.bam)?;
    bam.set_threads(8)?;

    for bed_rec in bed_records {
        bam.fetch((&bed_rec.chrom, bed_rec.position, bed_rec.position + 1))
            .unwrap_or_else(|_| {
                panic!(
                    "Failed to fetch region: {}:{}-{}",
                    &bed_rec.chrom,
                    bed_rec.position,
                    bed_rec.position + 1
                )
            });

        let records: Vec<bam::Record> = bam.records().map(|r| r.unwrap()).collect();
        let fibers = FiberseqData::from_records(records, &header_view, 0);
        for fiber in fibers {
            out.write(&fiber.record)?;
        }
    }

    Ok(())
}

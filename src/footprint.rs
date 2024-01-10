use crate::center::CenterPosition;

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
    modules: Vec<(i64, i64)>,
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

    pub fn max_pos(&self) -> i64 {
        self.modules[self.modules.len() - 1].1
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

use super::bam_writer;
use super::cli::FootprintOptions;
use super::fiber::*;
use anyhow;
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
}

pub fn start_finding_footprints(opts: &FootprintOptions) -> Result<(), anyhow::Error> {
    let yaml_buff = bio_io::buffer_from(&opts.yaml)?;
    let yaml: FootprintYaml = serde_yaml::from_reader(yaml_buff)?;
    yaml.check_for_valid_input()?;
    log::debug!("YAML: {:?}", yaml);

    let mut bam = bio_io::bam_reader(&opts.bam, 8);
    let mut out = bam_writer(&opts.out, &bam, 8);
    for mut rec in FiberseqRecords::new(&mut bam, 0) {
        // TODO
        rec.m6a.starts = vec![];
        out.write(&rec.record)?;
    }
    Ok(())
}

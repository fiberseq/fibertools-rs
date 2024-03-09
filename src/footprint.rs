use std::io::BufRead;

use super::bio_io;
use super::cli::FootprintOptions;
use super::fiber::*;
use crate::center::CenterPosition;
use anyhow;
//use rayon::prelude::*;
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
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: char,
    footprint: &'a FootprintYaml,
}

impl<'a> ReferenceMotif<'a> {
    pub fn new(line: &str, footprint: &'a FootprintYaml) -> anyhow::Result<Self> {
        let tokens = line.split('\t').collect::<Vec<_>>();
        assert!(tokens.len() >= 3);
        let st = tokens[1].parse::<i64>().unwrap();
        // subtract 1 from the end to make it inclusive ranges []
        let en = tokens[2].parse::<i64>().unwrap() - 1;
        anyhow::ensure!(
            st <= en,
            "BED start must be less than or equal to the end: {}",
            line
        );

        // get the strand for the 6 or 4th column
        let strand =
            if (tokens.len() >= 6 && tokens[5] == "-") || (tokens.len() >= 4 && tokens[3] == "-") {
                '-'
            } else {
                '+'
            };

        anyhow::ensure!(
            en - st + 1 == footprint.max_pos() as i64,
            "Motif length in the BED record ({}) does not match to the total length ({}) in the footprint YAML:\n{}",
            en - st + 1,
            footprint.max_pos(),
            line
        );

        Ok(Self {
            chrom: tokens[0].to_string(),
            start: st,
            end: en,
            strand,
            footprint,
        })
    }

    pub fn spans(&self, start: i64, end: i64) -> bool {
        // check if coordinates span the motif
        self.start > start && self.end < end
    }

    pub fn overlaps(&self, start: i64, end: i64) -> bool {
        // check if coordinates overlap the motif
        self.start < end && self.end > start
    }
}

pub struct Footprint<'a> {
    pub n_spanning_fibers: usize,
    pub n_spanning_msps: usize,
    pub has_spanning_msp: Vec<bool>,
    pub msp_fire_scores: Vec<i16>,
    pub n_overlapping_nucs: usize,
    pub has_overlapping_nucleosome: Vec<bool>,
    pub footprint_codes: Vec<u16>,
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
            msp_fire_scores: vec![],
            n_overlapping_nucs: 0,
            has_overlapping_nucleosome: vec![],
            footprint_codes: vec![],
            motif,
            fibers,
        };
        // add the number of msps that span the footprint
        footprint.spanning_msps();
        footprint.establish_footprints();
        footprint.overlapping_nucleosomes();
        footprint
    }

    fn spanning_msps(&mut self) {
        for fiber in self.fibers.iter() {
            let mut has_spanning_msp = false;
            let mut msp_qual = -1;
            for msp in &fiber.msp {
                // skip if there is no mapping of the msp
                match msp.4 {
                    Some((rs, re, _rl)) => {
                        if self.motif.spans(rs, re) {
                            self.n_spanning_msps += 1;
                            has_spanning_msp = true;
                            msp_qual = msp.3 as i16;
                            break;
                        }
                    }
                    None => continue,
                }
            }
            self.has_spanning_msp.push(has_spanning_msp);
            self.msp_fire_scores.push(msp_qual);
        }
    }

    pub fn overlapping_nucleosomes(&mut self) {
        // for each fiber, check if there is a nucleosome that overlaps with the motif
        for fiber in self.fibers.iter() {
            let mut has_overlapping_nucleosome = false;
            for nuc in &fiber.nuc {
                match nuc.4 {
                    Some((rs, re, _rl)) => {
                        if self.motif.overlaps(rs, re) {
                            self.n_overlapping_nucs += 1;
                            has_overlapping_nucleosome = true;
                            break;
                        }
                    }
                    None => continue,
                }
            }
            self.has_overlapping_nucleosome
                .push(has_overlapping_nucleosome);
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
        let mut m6a_vec = vec![false; motif_end];

        for m6a in fiber.m6a.reference_starts.iter().flatten() {
            if m6a < &self.motif.start {
                continue;
            }
            if m6a > &self.motif.end {
                break;
            }
            let m6a_rel_pos = m6a - self.motif.start;
            m6a_vec[m6a_rel_pos as usize] = true;
        }
        if self.motif.strand == '-' {
            m6a_vec.reverse();
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

    pub fn is_footprinted(&self, module: usize) -> Vec<bool> {
        let mut rtn = vec![];
        for footprint_code in &self.footprint_codes {
            if footprint_code & (1 << (module + 1)) != 0 {
                rtn.push(true);
            } else {
                rtn.push(false);
            }
        }
        rtn
    }

    pub fn footprinted_by_module_count(&self) -> Vec<usize> {
        let mut out = vec![];
        for i in 0..self.motif.footprint.modules.len() {
            let count = self.is_footprinted(i).iter().filter(|&&x| x).count();
            out.push(count);
        }
        out
    }

    pub fn out_bed_header(&self) -> String {
        let mut out =
            "#chrom\tstart\tend\tstrand\tn_spanning_fibers\tn_spanning_msps\tn_overlapping_nucs\t"
                .to_string();
        out += &self
            .motif
            .footprint
            .modules
            .iter()
            .map(|(st, en)| format!("module:{}-{}", st, en))
            .collect::<Vec<_>>()
            .join("\t");
        out += "\tfootprint_codes\tfire_quals\tfiber_names";
        out
    }
}

impl std::fmt::Display for Footprint<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut out = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.motif.chrom,
            self.motif.start,
            self.motif.end + 1,
            self.motif.strand,
            self.n_spanning_fibers,
            self.n_spanning_msps,
            self.n_overlapping_nucs,
        );

        out += &(self
            .footprinted_by_module_count()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join("\t"));

        out += &format!(
            "\t{}\t{}\t{}",
            self.footprint_codes
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
            self.msp_fire_scores
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
            self.fibers
                .iter()
                .map(|x| String::from_utf8_lossy(x.record.qname()))
                .collect::<Vec<_>>()
                .join(","),
        );
        writeln!(f, "{}", out)
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

    let bam = bio_io::bam_reader(&opts.bam, 1);
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let mut out = bio_io::writer(&opts.out)?;
    let mut bam = rust_htslib::bam::IndexedReader::from_path(&opts.bam)?;
    bam.set_threads(8)?;

    let reader = bio_io::buffer_from(&opts.bed)?;
    let mut first = true;
    // read lines in opts.bed
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let motif = ReferenceMotif::new(&line, &yaml)?;
        bam.fetch((&motif.chrom, motif.start, motif.end))?;
        let records: Vec<bam::Record> = bam.records().map(|r| r.unwrap()).collect();
        let fibers = FiberseqData::from_records(records, &header_view, 0);

        let footprint = Footprint::new(&motif, &fibers);
        if first {
            writeln!(out, "{}", footprint.out_bed_header())?;
            first = false;
        }
        out.write_all(format!("{}", footprint).as_bytes())?;
    }
    Ok(())
}

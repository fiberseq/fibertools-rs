use super::bamlift::*;
use super::extract::*;
use super::*;
use bio::alphabets::dna::revcomp;
use rust_htslib::bam::record;
use rust_htslib::bam::Read;
use rust_htslib::{bam, bam::ext::BamRecordExtensions};
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

#[derive(Clone)]
pub struct CenterPosition {
    pub chrom: String,
    pub position: i64,
    pub strand: char,
}
pub struct CenteredFiberData {
    fiber: FiberseqData,
    record: record::Record,
    offset: i64,
    center_position: CenterPosition,
}

impl CenteredFiberData {
    pub fn new(
        fiber: FiberseqData,
        record: bam::Record,
        center_position: CenterPosition,
    ) -> Option<Self> {
        let offset = CenteredFiberData::find_offset(&record, center_position.position);
        offset.map(|offset| CenteredFiberData {
            fiber,
            record,
            offset,
            center_position,
        })
    }

    pub fn get_centering_position(&self) -> i64 {
        self.center_position.position
    }

    /// find the query position that corresponds to the central reference position
    fn find_offset(record: &bam::Record, reference_position: i64) -> Option<i64> {
        let read_center = get_exact_query_positions(record, &[reference_position]);
        log::debug!(
            "{}, {}, {}, {:?}",
            reference_position,
            record.reference_start(),
            record.reference_end(),
            read_center
        );
        if read_center.is_empty() {
            None
        } else {
            Some(read_center[0])
        }
    }

    /// Center positions on the read around the reference position.
    fn apply_offset(&self, positions: &[i64]) -> Vec<i64> {
        let out = positions.iter().map(|&p| p - self.offset).collect();
        if self.center_position.strand == '+' {
            out
        } else {
            out.iter().rev().map(|&p| -p).collect()
        }
    }

    pub fn m6a_positions(&self) -> Vec<i64> {
        self.apply_offset(&self.fiber.base_mods.m6a_positions(false))
    }

    pub fn cpg_positions(&self) -> Vec<i64> {
        self.apply_offset(&self.fiber.base_mods.cpg_positions(false))
    }

    // todo
    pub fn nuc_positions(&self) -> Vec<(i64, i64)> {
        let starts = self.apply_offset(&self.fiber.get_nuc(false, true));
        let ends = self.apply_offset(
            &self
                .fiber
                .get_nuc(false, false)
                .iter()
                .zip(starts.iter())
                .map(|(length, start)| length + start)
                .collect::<Vec<_>>(),
        );
        starts
            .into_iter()
            .zip(ends)
            .map(|(s, e)| if s < e { (s, e) } else { (e, s) })
            .collect()
    }

    // todo
    pub fn msp_positions(&self) -> Vec<(i64, i64)> {
        let starts = self.apply_offset(&self.fiber.get_msp(false, true));
        let ends = self.apply_offset(
            &self
                .fiber
                .get_msp(false, false)
                .iter()
                .zip(starts.iter())
                .map(|(length, start)| length + start)
                .collect::<Vec<_>>(),
        );
        starts
            .into_iter()
            .zip(ends)
            .map(|(s, e)| if s < e { (s, e) } else { (e, s) })
            .collect()
    }

    pub fn get_sequence(&self) -> String {
        let forward_bases = if self.center_position.strand == '+' {
            self.record.seq().as_bytes()
        } else {
            revcomp(self.record.seq().as_bytes())
        };
        String::from_utf8_lossy(&forward_bases).to_string()
    }

    pub fn leading_columns(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.center_position.chrom,
            self.center_position.position,
            self.center_position.strand,
            self.record.reference_start(),
            self.record.reference_end(),
            std::str::from_utf8(self.record.qname()).unwrap(),
            -self.offset,
            self.record.seq_len() as i64 - self.offset,
            self.record.seq_len()
        )
    }

    // TODO
    pub fn write(&self) -> String {
        let m6a = join_by_str(self.m6a_positions(), ",");
        let cpg = join_by_str(self.cpg_positions(), ",");
        format!(
            "{}{}\t{}\t{}\n",
            self.leading_columns(),
            m6a,
            cpg,
            self.get_sequence(),
        )
    }
    pub fn header() -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "chrom",
            "centering_position",
            "strand",
            "reference_start",
            "reference_end",
            "query_name",
            "centered_query_start",
            "centered_query_end",
            "query_length",
            "centered_m6a_positions",
            "centered_cpg_positions",
            "query_sequence"
        )
    }
}

pub fn center(records: Vec<bam::Record>, center_position: CenterPosition) {
    let fiber_data = FiberseqData::from_records(&records);
    let iter = fiber_data.into_iter().zip(records.into_iter());
    let mut total = 0;
    let mut missing = 0;
    for (fiber, record) in iter {
        let out = match CenteredFiberData::new(fiber, record, center_position.clone()) {
            Some(centered_fiber) => {
                total += 1;
                centered_fiber.write()
            }
            None => {
                total += 1;
                missing += 1;
                "".to_string()
            }
        };
        print!("{}", out);
    }
    if missing > 0 {
        log::info!(
            "Unable to exactly map {}/{} reads at position {}:{}",
            missing,
            total,
            center_position.chrom.clone(),
            center_position.position
        );
    }
}

pub fn center_fiberdata(bam: &mut bam::IndexedReader, center_positions: Vec<CenterPosition>) {
    print!("{}", CenteredFiberData::header());
    for center_position in center_positions {
        bam.fetch((
            &center_position.chrom,
            center_position.position,
            center_position.position + 1,
        ))
        .expect("Failed to fetch region");
        let records: Vec<bam::Record> = bam.records().map(|r| r.unwrap()).collect();
        center(records, center_position);
    }
}

pub fn read_center_positions(infile: &str) -> io::Result<Vec<CenterPosition>> {
    let file = File::open(infile)?;
    let reader = BufReader::new(file);
    let mut rtn = vec![];
    for line in reader.lines() {
        let line = line?;
        let tokens = line.split('\t').collect::<Vec<_>>();
        assert!(tokens.len() >= 2);
        let strand = if tokens.len() >= 3 && tokens[2] == "-" {
            '-'
        } else {
            '+'
        };
        rtn.push(CenterPosition {
            chrom: tokens[0].to_string(),
            position: tokens[1].parse::<i64>().unwrap(),
            strand,
        });
    }
    Ok(rtn)
}

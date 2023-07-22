use super::extract::*;
use super::*;
use bamlift::*;
use bio::alphabets::dna::revcomp;
use bio_io;
use indicatif::{style, ProgressBar};
use rayon::prelude::*;
use rust_htslib::bam::Read;
use rust_htslib::{bam, bam::ext::BamRecordExtensions};
use std::fmt::Write;
use std::io::{self, prelude::*};

#[derive(Clone)]
pub struct CenterPosition {
    pub chrom: String,
    pub position: i64,
    pub strand: char,
}
pub struct CenteredFiberData {
    fiber: FiberseqData,
    rg: String,
    offset: i64,
    pub dist: Option<i64>,
    center_position: CenterPosition,
    pub reference: bool,
}

impl CenteredFiberData {
    pub fn new(
        fiber: FiberseqData,
        center_position: CenterPosition,
        dist: Option<i64>,
        reference: bool,
    ) -> Option<Self> {
        let offset = if reference {
            Some(center_position.position)
        } else {
            CenteredFiberData::find_offset(&fiber.record, center_position.position)
        };
        // get RG
        let rg = if let Ok(bam::record::Aux::String(f)) = fiber.record.aux(b"RG") {
            f
        } else {
            "."
        }
        .to_string();

        offset.map(|offset| CenteredFiberData {
            fiber,
            rg,
            offset,
            dist,
            center_position,
            reference,
        })
    }

    pub fn get_centering_position(&self) -> i64 {
        self.center_position.position
    }

    /// find the query position that corresponds to the central reference position
    fn find_offset(record: &bam::Record, reference_position: i64) -> Option<i64> {
        let read_center: Vec<i64> = lift_query_positions_exact(record, &[reference_position])
            .into_iter()
            .filter(|&x| x >= 0)
            .collect();
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

    /// Get the sequence
    pub fn subset_sequence(&self) -> String {
        let dist = if let Some(dist) = self.dist { dist } else { 0 };
        let seq = self.fiber.record.seq().as_bytes();

        let mut out_seq: Vec<u8> = vec![];
        let st = self.offset - dist; //(self.offset - dist)
        for pos in st..(self.offset + dist + 1) {
            if pos < 0 || pos as usize >= seq.len() {
                out_seq.push(b'N');
            } else {
                out_seq.push(seq[pos as usize]);
            }
        }
        if self.center_position.strand == '-' {
            out_seq = revcomp(out_seq);
        }
        //assert_eq!(out_seq.len() as i64, dist * 2 + 1);
        String::from_utf8_lossy(&out_seq).to_string()
    }

    pub fn m6a_positions(&self) -> Vec<i64> {
        self.apply_offset(&self.fiber.base_mods.m6a_positions(self.reference))
    }

    pub fn cpg_positions(&self) -> Vec<i64> {
        self.apply_offset(&self.fiber.base_mods.cpg_positions(self.reference))
    }

    fn get_start_end_positions(&self, starts: &[i64], ends: &[i64]) -> Vec<(i64, i64)> {
        let starts = self.apply_offset(starts);
        let ends = self.apply_offset(ends);
        starts
            .into_iter()
            .zip(ends)
            .map(|(s, e)| {
                if e - s > 300 || s - e > 300 {
                    log::trace!("{}, {}, {}", s, e, e - s);
                }
                if s < e {
                    (s, e)
                } else {
                    (e, s)
                }
            })
            .collect()
    }

    pub fn nuc_positions(&self) -> Vec<(i64, i64)> {
        let (starts, ends) = if self.reference {
            (
                &self.fiber.nuc.reference_starts,
                &self.fiber.nuc.reference_ends,
            )
        } else {
            (&self.fiber.nuc.starts, &self.fiber.nuc.ends)
        };
        self.get_start_end_positions(starts, ends)
    }

    pub fn msp_positions(&self) -> Vec<(i64, i64)> {
        let (starts, ends) = if self.reference {
            (
                &self.fiber.msp.reference_starts,
                &self.fiber.msp.reference_ends,
            )
        } else {
            (&self.fiber.msp.starts, &self.fiber.msp.ends)
        };
        self.get_start_end_positions(starts, ends)
    }

    pub fn get_sequence(&self) -> String {
        let forward_bases = if self.center_position.strand == '+' {
            self.fiber.record.seq().as_bytes()
        } else {
            revcomp(self.fiber.record.seq().as_bytes())
        };
        String::from_utf8_lossy(&forward_bases).to_string()
    }

    pub fn leading_columns(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.center_position.chrom,
            self.center_position.position,
            self.center_position.strand,
            self.subset_sequence(),
            self.fiber.record.reference_start(),
            self.fiber.record.reference_end(),
            std::str::from_utf8(self.fiber.record.qname()).unwrap(),
            self.rg,
            -self.offset,
            self.fiber.record.seq_len() as i64 - self.offset,
            self.fiber.record.seq_len()
        )
    }

    pub fn write(&self) -> String {
        let m6a = join_by_str(self.m6a_positions(), ",");
        let cpg = join_by_str(self.cpg_positions(), ",");
        let (nuc_st, nuc_en) = unzip_to_vectors(self.nuc_positions());
        let (msp_st, msp_en) = unzip_to_vectors(self.msp_positions());
        format!(
            "{}{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.leading_columns(),
            m6a,
            cpg,
            join_by_str(nuc_st, ","),
            join_by_str(nuc_en, ","),
            join_by_str(msp_st, ","),
            join_by_str(msp_en, ","),
            self.get_sequence(),
        )
    }

    pub fn leading_header() -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            "chrom",
            "centering_position",
            "strand",
            "subset_sequence",
            "reference_start",
            "reference_end",
            "query_name",
            "RG",
            "centered_query_start",
            "centered_query_end",
            "query_length",
        )
    }
    pub fn header() -> String {
        format!(
            "{}{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            CenteredFiberData::leading_header(),
            "centered_m6a_positions",
            "centered_5mC_positions",
            "centered_nuc_starts",
            "centered_nuc_ends",
            "centered_msp_starts",
            "centered_msp_ends",
            "query_sequence"
        )
    }

    pub fn long_header() -> String {
        format!(
            "{}{}\t{}\t{}\n",
            CenteredFiberData::leading_header(),
            "centered_position_type",
            "centered_start",
            "centered_end",
        )
    }

    pub fn write_long(&self) -> String {
        let m6a = self.m6a_positions();
        let cpg = self.cpg_positions();
        let (nuc_st, nuc_en) = unzip_to_vectors(self.nuc_positions());
        let (msp_st, msp_en) = unzip_to_vectors(self.msp_positions());
        let mut rtn = String::new();
        for (t, vals) in [
            ("m6a", (m6a, None)),
            ("5mC", (cpg, None)),
            ("nuc", (nuc_st, Some(nuc_en))),
            ("msp", (msp_st, Some(msp_en))),
        ] {
            let starts = vals.0;
            let ends = match vals.1 {
                Some(ends) => ends,
                None => starts.iter().map(|&st| st + 1).collect(),
            };
            for (&st, &en) in starts.iter().zip(ends.iter()) {
                let mut write = true;
                if let Some(dist) = self.dist {
                    // skip writing if we are outside the motif range
                    if en <= -dist || st > dist {
                        write = false;
                    }
                };
                if write {
                    // add the leading data
                    rtn.push_str(&self.leading_columns());
                    // add the long form data
                    write!(&mut rtn, "{}", format_args!("{}\t{}\t{}\n", t, st, en)).unwrap();
                }
            }
        }

        rtn
    }
}

#[allow(clippy::too_many_arguments)]
pub fn center(
    records: Vec<bam::Record>,
    header_view: &rust_htslib::bam::HeaderView,
    center_position: CenterPosition,
    min_ml_score: u8,
    wide: bool,
    dist: Option<i64>,
    reference: bool,
    buffer: &mut Box<dyn std::io::Write>,
) {
    let fiber_data = FiberseqData::from_records(records, header_view, min_ml_score);
    let total = fiber_data.len();
    let mut seen = 0;

    fiber_data
        .into_par_iter()
        .map(|fiber| {
            match CenteredFiberData::new(fiber, center_position.clone(), dist, reference) {
                Some(centered_fiber) => {
                    if wide {
                        centered_fiber.write()
                    } else {
                        centered_fiber.write_long()
                    }
                }
                None => "".to_string(),
            }
        })
        .collect::<Vec<_>>()
        .iter()
        .filter(|x| !x.is_empty())
        .for_each(|x| {
            seen += 1;
            write_to_file(x, buffer);
        });

    if total - seen > 1 {
        log::warn!(
            "Unable to exactly map {}/{} reads at position {}:{}",
            total - seen,
            total,
            center_position.chrom.clone(),
            center_position.position
        );
    }
}

pub fn center_fiberdata(
    bam: &mut bam::IndexedReader,
    center_positions: Vec<CenterPosition>,
    min_ml_score: u8,
    wide: bool,
    dist: Option<i64>,
    reference: bool,
) {
    // header needed for the contig name...
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    //
    let mut buffer = bio_io::writer("-").unwrap();

    if wide {
        bio_io::write_to_file(&CenteredFiberData::header(), &mut buffer);
    } else {
        bio_io::write_to_file(&CenteredFiberData::long_header(), &mut buffer);
    }

    let pb = ProgressBar::new(center_positions.len() as u64);
    pb.set_style(
        style::ProgressStyle::with_template(PROGRESS_STYLE)
            .unwrap()
            .progress_chars("##-"),
    );

    for center_position in center_positions {
        bam.fetch((
            &center_position.chrom,
            center_position.position,
            center_position.position + 1,
        ))
        .expect("Failed to fetch region");
        let records: Vec<bam::Record> = bam.records().map(|r| r.unwrap()).collect();
        center(
            records,
            &header_view,
            center_position,
            min_ml_score,
            wide,
            dist,
            reference,
            &mut buffer,
        );
        pb.inc(1);
    }
    buffer.flush().unwrap();
    pb.finish_with_message("\ndone");
}

pub fn read_center_positions(infile: &str) -> io::Result<Vec<CenterPosition>> {
    let reader = bio_io::buffer_from(infile).expect("Failed to open bed file");
    let mut rtn = vec![];
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let tokens = line.split('\t').collect::<Vec<_>>();
        assert!(tokens.len() >= 3);
        let st = tokens[1].parse::<i64>().unwrap();
        let en = tokens[2].parse::<i64>().unwrap();
        // get the strand for the 6 or 4th column
        let strand =
            if (tokens.len() >= 6 && tokens[5] == "-") || (tokens.len() >= 4 && tokens[3] == "-") {
                '-'
            } else {
                '+'
            };

        let (strand, position) = if strand == '-' {
            ('-', en - 1)
        } else {
            ('+', st)
        };
        rtn.push(CenterPosition {
            chrom: tokens[0].to_string(),
            position,
            strand,
        });
    }
    Ok(rtn)
}

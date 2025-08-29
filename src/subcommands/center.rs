use crate::cli::CenterOptions;
use crate::fiber::FiberseqData;
use crate::utils::bamlift::*;
use crate::utils::bio_io;
use crate::*;
use bio::alphabets::dna::revcomp;
use indicatif::{style, ProgressBar};
use rayon::prelude::*;
use rust_htslib::bam::Read;
use rust_htslib::{bam, bam::ext::BamRecordExtensions};
use std::fmt::Write;
use std::io::{self, prelude::*};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CenterPosition {
    pub chrom: String,
    pub position: i64,
    pub strand: char,
}
pub struct CenteredFiberData {
    fiber: FiberseqData,
    pub dist: Option<i64>,
    center_position: CenterPosition,
    pub offset: i64,
    pub reference: bool,
    pub simplify: bool,
}

impl CenteredFiberData {
    pub fn new(
        fiber: FiberseqData,
        center_position: CenterPosition,
        dist: Option<i64>,
        reference: bool,
        simplify: bool,
    ) -> Option<Self> {
        let (ref_offset, mol_offset) =
            CenteredFiberData::find_offsets(&fiber.record, &center_position);
        let offset = if reference { ref_offset } else { mol_offset };

        let fiber = fiber.center(&center_position)?;

        Some(CenteredFiberData {
            fiber,
            dist,
            center_position,
            offset,
            reference,
            simplify,
        })
    }
    /// find both the ref and mol offsets
    pub fn find_offsets(record: &bam::Record, center_position: &CenterPosition) -> (i64, i64) {
        let ref_offset = center_position.position;
        let mol_offset =
            CenteredFiberData::find_offset(record, center_position.position).unwrap_or(0);
        (ref_offset, mol_offset)
    }

    /// find the query position that corresponds to the central reference position
    pub fn find_offset(record: &bam::Record, reference_position: i64) -> Option<i64> {
        let read_center: Vec<i64> = lift_query_positions_exact(record, &[reference_position])
            .ok()?
            .into_iter()
            .flatten()
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

    /// Get the sequence
    pub fn subset_sequence(&self) -> String {
        if self.simplify {
            return "N".to_string();
        }
        let dist = self.dist.unwrap_or(0);
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

    pub fn get_sequence(&self) -> String {
        let forward_bases = if self.center_position.strand == '+' {
            self.fiber.record.seq().as_bytes()
        } else {
            revcomp(self.fiber.record.seq().as_bytes())
        };
        String::from_utf8_lossy(&forward_bases).to_string()
    }

    pub fn leading_columns(&self) -> String {
        let (mut c_query_start, mut c_query_end) = if self.reference {
            (
                self.fiber.record.reference_start() - self.center_position.position,
                self.fiber.record.reference_end() - self.center_position.position,
            )
        } else {
            let query_length = self.fiber.record.seq_len() as i64;
            (-self.offset, query_length - self.offset)
        };

        if self.center_position.strand == '-' {
            c_query_start = -c_query_start;
            c_query_end = -c_query_end;
            if c_query_start > c_query_end {
                std::mem::swap(&mut c_query_start, &mut c_query_end);
            }
        }
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.center_position.chrom,
            self.center_position.position,
            self.center_position.strand,
            self.subset_sequence(),
            self.fiber.record.reference_start(),
            self.fiber.record.reference_end(),
            std::str::from_utf8(self.fiber.record.qname()).unwrap(),
            self.fiber.rg,
            self.fiber.get_hp(),
            c_query_start,
            c_query_end,
            self.fiber.record.seq_len()
        )
    }

    #[allow(clippy::type_complexity)]
    fn grab_data(
        &self,
    ) -> (
        Vec<Option<i64>>,
        Vec<u8>,
        Vec<Option<i64>>,
        Vec<u8>,
        Vec<Option<i64>>,
        Vec<Option<i64>>,
        Vec<Option<i64>>,
        Vec<Option<i64>>,
        Vec<u8>,
    ) {
        if self.reference {
            (
                self.fiber.m6a.reference_starts(),
                self.fiber.m6a.qual(),
                self.fiber.cpg.reference_starts(),
                self.fiber.cpg.qual(),
                self.fiber.nuc.reference_starts(),
                self.fiber.nuc.reference_ends(),
                self.fiber.msp.reference_starts(),
                self.fiber.msp.reference_ends(),
                self.fiber.msp.qual(),
            )
        } else {
            (
                self.fiber.m6a.option_starts(),
                self.fiber.m6a.qual(),
                self.fiber.cpg.option_starts(),
                self.fiber.cpg.qual(),
                self.fiber.nuc.option_starts(),
                self.fiber.nuc.option_ends(),
                self.fiber.msp.option_starts(),
                self.fiber.msp.option_ends(),
                self.fiber.msp.qual(),
            )
        }
    }

    pub fn write(&self) -> String {
        let (m6a, m6a_qual, cpg, cpg_qual, nuc_st, nuc_en, msp_st, msp_en, fire) = self.grab_data();
        format!(
            "{}{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.leading_columns(),
            join_by_str_option(&m6a, ","),
            join_by_str(&m6a_qual, ","),
            join_by_str_option(&cpg, ","),
            join_by_str(&cpg_qual, ","),
            join_by_str_option(&nuc_st, ","),
            join_by_str_option(&nuc_en, ","),
            join_by_str_option(&msp_st, ","),
            join_by_str_option(&msp_en, ","),
            join_by_str(&fire, ","),
            self.get_sequence(),
        )
    }

    pub fn leading_header() -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            "chrom",
            "centering_position",
            "strand",
            "subset_sequence",
            "reference_start",
            "reference_end",
            "query_name",
            "RG",
            "HP",
            "centered_query_start",
            "centered_query_end",
            "query_length",
        )
    }
    pub fn header() -> String {
        format!(
            "{}{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            CenteredFiberData::leading_header(),
            "centered_m6a_positions",
            "m6a_qual",
            "centered_5mC_positions",
            "5mC_qual",
            "centered_nuc_starts",
            "centered_nuc_ends",
            "centered_msp_starts",
            "centered_msp_ends",
            "query_sequence"
        )
    }

    pub fn long_header() -> String {
        format!(
            "{}{}\t{}\t{}\t{}\n",
            CenteredFiberData::leading_header(),
            "centered_position_type",
            "centered_start",
            "centered_end",
            "centered_qual",
        )
    }

    pub fn write_long(&self) -> String {
        let mut rtn = String::new();
        let (m6a, m6a_qual, cpg, cpg_qual, nuc_st, nuc_en, msp_st, msp_en, fire) = self.grab_data();
        for (t, vals) in [
            ("m6a", (m6a, None, Some(m6a_qual))),
            ("5mC", (cpg, None, Some(cpg_qual))),
            ("nuc", (nuc_st, Some(nuc_en), None)),
            ("msp", (msp_st, Some(msp_en), Some(fire))),
        ] {
            let starts = vals.0.iter().collect::<Vec<_>>();
            let ends: Vec<Option<i64>> = match vals.1 {
                Some(ends) => ends.to_vec(),
                None => vec![None; starts.len()],
            };
            let quals = match vals.2 {
                Some(quals) => quals.to_vec(),
                None => vec![0; starts.len()],
            };
            let mut write_count = 0;
            for ((&st, &en), &qual) in starts.iter().zip(ends.iter()).zip(quals.iter()) {
                let Some(st) = st else {
                    continue;
                };
                let st = *st;
                let en = en.unwrap_or(st + 1);

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
                    write!(
                        &mut rtn,
                        "{}",
                        format_args!("{}\t{}\t{}\t{}\n", t, st, en, qual)
                    )
                    .unwrap();
                    write_count += 1;
                }
            }
            log::debug!("{t}: {write_count}");
        }

        rtn
    }
}

#[allow(clippy::too_many_arguments)]
pub fn center(
    records: Vec<bam::Record>,
    header_view: &rust_htslib::bam::HeaderView,
    center_position: CenterPosition,
    opts: &CenterOptions,
    buffer: &mut Box<dyn std::io::Write>,
) {
    let fiber_data = FiberseqData::from_records(records, header_view, &opts.input.filters);
    let total = fiber_data.len();
    let mut seen = 0;

    let to_write: Vec<String> = fiber_data
        .into_par_iter()
        .map(|fiber| {
            match CenteredFiberData::new(
                fiber,
                center_position.clone(),
                opts.dist,
                opts.reference,
                opts.simplify,
            ) {
                Some(centered_fiber) => {
                    if opts.wide {
                        centered_fiber.write()
                    } else {
                        centered_fiber.write_long()
                    }
                }
                None => "".to_string(),
            }
        })
        .filter(|x| !x.is_empty())
        .collect::<Vec<_>>();

    for line in to_write {
        seen += 1;
        write_to_file(&line, buffer);
    }

    log::debug!(
        "centering {} records of {} on {}:{}:{}",
        seen,
        total,
        center_position.chrom,
        center_position.position,
        center_position.strand
    );

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

pub fn center_fiberdata(center_opts: &mut CenterOptions) -> anyhow::Result<()> {
    let mut bam = center_opts.input.indexed_bam_reader();
    let center_positions = read_center_positions(&center_opts.bed)?;

    // header needed for the contig name...
    let header_view = center_opts.input.header_view();
    // output buffer
    let mut buffer = bio_io::writer("-").unwrap();

    if center_opts.wide {
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
        .unwrap_or_else(|_| {
            panic!(
                "Failed to fetch region: {}:{}-{}",
                &center_position.chrom,
                center_position.position,
                center_position.position + 1
            )
        });

        let records: Vec<bam::Record> = center_opts
            .input
            .filters
            .filter_on_bit_flags(bam.records())
            .collect();

        center(
            records,
            &header_view,
            center_position,
            center_opts,
            &mut buffer,
        );
        pb.inc(1);
    }
    buffer.flush().unwrap();
    pb.finish_with_message("\ndone");
    Ok(())
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

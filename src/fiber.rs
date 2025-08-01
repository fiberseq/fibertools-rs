use super::subcommands::center::CenterPosition;
use super::subcommands::center::CenteredFiberData;
use super::utils::input_bam::FiberFilters;
use super::*;
use crate::utils::bamannotations::*;
use crate::utils::basemods::BaseMods;
use crate::utils::bio_io::*;
use crate::utils::ftexpression::apply_filter_fsd;
use rayon::prelude::*;
use rust_htslib::bam::Read;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux, bam::HeaderView};
use std::collections::HashMap;
use std::fmt::Write;

#[derive(Debug, Clone, PartialEq)]
pub struct FiberseqData {
    pub record: bam::Record,
    pub msp: Ranges,
    pub nuc: Ranges,
    pub m6a: Ranges,
    pub cpg: Ranges,
    pub base_mods: BaseMods,
    pub ec: f32,
    pub target_name: String,
    pub rg: String,
    pub center_position: Option<CenterPosition>,
}

impl FiberseqData {
    pub fn new(record: bam::Record, target_name: Option<&String>, filters: &FiberFilters) -> Self {
        // read group
        let rg = if let Ok(Aux::String(f)) = record.aux(b"RG") {
            log::trace!("{f}");
            f
        } else {
            "."
        }
        .to_string();

        let nuc_starts = get_u32_tag(&record, b"ns");
        let msp_starts = get_u32_tag(&record, b"as");
        let nuc_length = get_u32_tag(&record, b"nl");
        let msp_length = get_u32_tag(&record, b"al");
        let nuc = Ranges::new(&record, nuc_starts, None, Some(nuc_length));
        let mut msp = Ranges::new(&record, msp_starts, None, Some(msp_length));
        let msp_qual = get_u8_tag(&record, b"aq");
        if !msp_qual.is_empty() {
            msp.set_qual(msp_qual);
        }

        // get the number of passes
        let ec = if let Ok(Aux::Float(f)) = record.aux(b"ec") {
            log::trace!("{f}");
            f
        } else {
            0.0
        };

        let target_name = match target_name {
            Some(t) => t.clone(),
            None => ".".to_string(),
        };

        // get fiberseq basemods
        let mut base_mods = BaseMods::new(&record, filters.min_ml_score);
        base_mods.filter_at_read_ends(filters.strip_starting_basemods);

        //let (m6a, cpg) = FiberMods::new(&base_mods);
        let m6a = base_mods.m6a();
        let cpg = base_mods.cpg();

        let mut fsd = FiberseqData {
            record,
            msp,
            nuc,
            m6a,
            base_mods,
            cpg,
            ec,
            target_name,
            rg,
            center_position: None,
        };

        apply_filter_fsd(&mut fsd, filters).expect("Failed to apply filter to FiberseqData");
        fsd
    }

    pub fn dict_from_head_view(head_view: &HeaderView) -> HashMap<i32, String> {
        if head_view.target_count() == 0 {
            return HashMap::new();
        }
        let target_u8s = head_view.target_names();
        let tids = target_u8s
            .iter()
            .map(|t| head_view.tid(t).expect("Unable to get tid"));
        let target_names = target_u8s
            .iter()
            .map(|&a| String::from_utf8_lossy(a).to_string());

        tids.zip(target_names)
            .map(|(id, t)| (id as i32, t))
            .collect()
    }

    pub fn target_name_from_tid(tid: i32, target_dict: &HashMap<i32, String>) -> Option<&String> {
        target_dict.get(&tid)
    }

    pub fn from_records(
        records: Vec<bam::Record>,
        head_view: &HeaderView,
        filters: &FiberFilters,
    ) -> Vec<Self> {
        let target_dict = Self::dict_from_head_view(head_view);
        records
            .into_par_iter()
            .map(|r| {
                let tid = r.tid();
                (r, Self::target_name_from_tid(tid, &target_dict))
            })
            .map(|(r, target_name)| Self::new(r, target_name, filters))
            .collect::<Vec<_>>()
    }

    //
    // GET FUNCTIONS
    //

    pub fn get_qname(&self) -> String {
        String::from_utf8_lossy(self.record.qname()).to_string()
    }

    pub fn get_rq(&self) -> Option<f32> {
        if let Ok(Aux::Float(f)) = self.record.aux(b"rq") {
            Some(f)
        } else {
            None
        }
    }

    pub fn get_hp(&self) -> String {
        if let Ok(Aux::U8(f)) = self.record.aux(b"HP") {
            format!("H{f}")
        } else {
            "UNK".to_string()
        }
    }

    /// Center all coordinates on the read using the offset attribute.
    pub fn center(&self, center_position: &CenterPosition) -> Option<Self> {
        // setup new fiberseq data object to return
        let mut new = self.clone();
        let (ref_offset, mol_offset) =
            CenteredFiberData::find_offsets(&self.record, center_position);

        // Apply offsets to all annotations using the new methods
        new.m6a
            .apply_offset(mol_offset, ref_offset, center_position.strand);
        new.cpg
            .apply_offset(mol_offset, ref_offset, center_position.strand);
        new.msp
            .apply_offset(mol_offset, ref_offset, center_position.strand);
        new.nuc
            .apply_offset(mol_offset, ref_offset, center_position.strand);

        // Validate that MSPs still start and end on m6A marks after centering
        new.validate_msp_m6a_alignment();

        Some(new)
    }

    /// Validate that all MSP boundaries align with m6A positions after centering
    fn validate_msp_m6a_alignment(&self) {
        let m6a_positions = self.m6a.starts();
        let msp_boundaries: Vec<i64> = self
            .msp
            .starts()
            .into_iter()
            .chain(self.msp.ends().into_iter().map(|x| x - 1))
            .collect();

        if m6a_positions.is_empty() || msp_boundaries.is_empty() {
            return; // Skip validation if no data
        }

        for msp_pos in &msp_boundaries {
            if !m6a_positions.contains(msp_pos) {
                log::warn!(
                    "MSP boundary at position {} does not align with m6A mark after centering in read {}",
                    msp_pos,
                    String::from_utf8_lossy(self.record.qname())
                );
            }
        }
    }

    //
    //  WRITE BED12 FUNCTIONS
    //
    pub fn write_msp(&self, reference: bool) -> String {
        let (starts, _ends, lengths) = if reference {
            (
                self.msp.reference_starts(),
                self.msp.reference_ends(),
                self.msp.reference_lengths(),
            )
        } else {
            (
                self.msp.option_starts(),
                self.msp.option_ends(),
                self.msp.option_lengths(),
            )
        };
        self.to_bed12(reference, &starts, &lengths, LINKER_COLOR)
    }

    pub fn write_nuc(&self, reference: bool) -> String {
        let (starts, _ends, lengths) = if reference {
            (
                self.nuc.reference_starts(),
                self.nuc.reference_ends(),
                self.nuc.reference_lengths(),
            )
        } else {
            (
                self.nuc.option_starts(),
                self.nuc.option_ends(),
                self.nuc.option_lengths(),
            )
        };
        self.to_bed12(reference, &starts, &lengths, NUC_COLOR)
    }

    pub fn write_m6a(&self, reference: bool) -> String {
        let starts = if reference {
            self.m6a.reference_starts()
        } else {
            self.m6a.option_starts()
        };
        let lengths = vec![Some(1); starts.len()];
        self.to_bed12(reference, &starts, &lengths, M6A_COLOR)
    }

    pub fn write_cpg(&self, reference: bool) -> String {
        let starts = if reference {
            self.cpg.reference_starts()
        } else {
            self.cpg.option_starts()
        };
        let lengths = vec![Some(1); starts.len()];
        self.to_bed12(reference, &starts, &lengths, CPG_COLOR)
    }

    pub fn to_bed12(
        &self,
        reference: bool,
        starts: &[Option<i64>],
        lengths: &[Option<i64>],
        color: &str,
    ) -> String {
        if starts.is_empty() {
            return "".to_string();
        }
        // skip if no alignments are here
        if self.record.is_unmapped() && reference {
            return "".to_string();
        }

        let ct;
        let start;
        let end;
        let name = String::from_utf8_lossy(self.record.qname()).to_string();
        let mut rtn: String = String::with_capacity(0);
        if reference {
            ct = &self.target_name;
            start = self.record.reference_start();
            end = self.record.reference_end();
        } else {
            ct = &name;
            start = 0;
            end = self.record.seq_len() as i64;
        }
        let score = self.ec.round() as i64;
        let strand = if self.record.is_reverse() { '-' } else { '+' };
        // filter out positions that do not have an exact liftover
        let (filtered_starts, filtered_lengths): (Vec<i64>, Vec<i64>) = starts
            .iter()
            .flatten()
            .zip(lengths.iter().flatten())
            .unzip();
        // skip empty ones
        if filtered_lengths.is_empty() || filtered_starts.is_empty() {
            return "".to_string();
        }
        let b_ct = filtered_starts.len() + 2;
        let b_ln: String = filtered_lengths
            .iter()
            .map(|&ln| ln.to_string() + ",")
            .collect();
        let b_st: String = filtered_starts
            .iter()
            .map(|&st| (st - start).to_string() + ",")
            .collect();
        assert_eq!(filtered_lengths.len(), filtered_starts.len());

        rtn.push_str(ct);
        rtn.push('\t');
        rtn.push_str(&start.to_string());
        rtn.push('\t');
        rtn.push_str(&end.to_string());
        rtn.push('\t');
        rtn.push_str(&name);
        rtn.push('\t');
        rtn.push_str(&score.to_string());
        rtn.push('\t');
        rtn.push(strand);
        rtn.push('\t');
        rtn.push_str(&start.to_string());
        rtn.push('\t');
        rtn.push_str(&end.to_string());
        rtn.push('\t');
        rtn.push_str(color);
        rtn.push('\t');
        rtn.push_str(&b_ct.to_string());
        rtn.push_str("\t0,"); // add a zero length start
        rtn.push_str(&b_ln);
        rtn.push_str("1\t0,"); // add a 1 base length and a 0 start point
        rtn.push_str(&b_st);
        write!(&mut rtn, "{}", format_args!("{}\n", end - start - 1)).unwrap();
        rtn
    }

    //
    // WRITE ALL FUNCTIONS
    //

    pub fn all_header(simplify: bool, quality: bool) -> String {
        let mut x = format!(
            "#{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            "ct", "st", "en", "fiber", "score", "strand", "sam_flag", "HP", "RG", "fiber_length",
        );
        if !simplify {
            x.push_str("fiber_sequence\t")
        }
        if quality {
            x.push_str("fiber_qual\t")
        }
        x.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "ec",
            "rq",
            "total_AT_bp",
            "total_m6a_bp",
            "total_nuc_bp",
            "total_msp_bp",
            "total_5mC_bp",
            "nuc_starts",
            "nuc_lengths",
            "ref_nuc_starts",
            "ref_nuc_lengths",
            "msp_starts",
            "msp_lengths",
            "fire",
            "ref_msp_starts",
            "ref_msp_lengths",
            "m6a",
            "ref_m6a",
            "m6a_qual",
            "5mC",
            "ref_5mC",
            "5mC_qual"
        ));
        x
    }

    pub fn write_all(&self, simplify: bool, quality: bool) -> String {
        // PB features
        let name = std::str::from_utf8(self.record.qname()).unwrap();
        let score = self.ec.round() as i64;
        let q_len = self.record.seq_len() as i64;
        let rq = match self.get_rq() {
            Some(x) => format!("{x}"),
            None => ".".to_string(),
        };
        // reference features
        let ct;
        let start;
        let end;
        let strand;
        if self.record.is_unmapped() {
            ct = ".";
            start = 0;
            end = 0;
            strand = '.';
        } else {
            ct = &self.target_name;
            start = self.record.reference_start();
            end = self.record.reference_end();
            strand = if self.record.is_reverse() { '-' } else { '+' };
        }
        let sam_flag = self.record.flags();
        let hp = self.get_hp();

        let at_count = self
            .record
            .seq()
            .as_bytes()
            .iter()
            .filter(|&x| *x == b'A' || *x == b'T')
            .count() as i64;

        // get the info
        let m6a_count = self.m6a.annotations.len();
        let m6a_qual = self.m6a.qual().iter().map(|a| Some(*a as i64)).collect();
        let cpg_count = self.cpg.annotations.len();
        let cpg_qual = self.cpg.qual().iter().map(|a| Some(*a as i64)).collect();
        let fire = self.msp.qual().iter().map(|a| Some(*a as i64)).collect();

        // write the features
        let mut rtn = String::with_capacity(0);
        // add first things 7
        rtn.write_fmt(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            ct, start, end, name, score, strand, sam_flag, hp, self.rg, q_len
        ))
        .unwrap();
        // add sequence
        if !simplify {
            rtn.write_fmt(format_args!(
                "{}\t",
                String::from_utf8_lossy(&self.record.seq().as_bytes()),
            ))
            .unwrap();
        }
        if quality {
            // TODO add quality offset
            rtn.write_fmt(format_args!(
                "{}\t",
                String::from_utf8_lossy(
                    &self
                        .record
                        .qual()
                        .iter()
                        .map(|x| x + 33)
                        .collect::<Vec<u8>>()
                ),
            ))
            .unwrap();
        }
        // add PB features
        let total_nuc_bp = self.nuc.lengths().iter().sum::<i64>();
        let total_msp_bp = self.msp.lengths().iter().sum::<i64>();
        rtn.write_fmt(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.ec, rq, at_count, m6a_count, total_nuc_bp, total_msp_bp, cpg_count
        ))
        .unwrap();
        // add fiber features
        let vecs = [
            self.nuc.option_starts(),
            self.nuc.option_lengths(),
            self.nuc.reference_starts(),
            self.nuc.reference_lengths(),
            self.msp.option_starts(),
            self.msp.option_lengths(),
            fire,
            self.msp.reference_starts(),
            self.msp.reference_lengths(),
            self.m6a.option_starts(),
            self.m6a.reference_starts(),
            m6a_qual,
            self.cpg.option_starts(),
            self.cpg.reference_starts(),
            cpg_qual,
        ];
        for vec in &vecs {
            if vec.is_empty() {
                rtn.push('.');
                rtn.push('\t');
            } else {
                let z: String = vec
                    .iter()
                    .map(|x| match x {
                        Some(y) => *y,
                        None => -1,
                    })
                    .map(|x| x.to_string() + ",")
                    .collect();
                rtn.write_fmt(format_args!("{z}\t")).unwrap();
            }
        }
        // replace the last tab with a newline
        let len = rtn.len();
        rtn.replace_range(len - 1..len, "\n");

        rtn
    }
}

pub struct FiberseqRecords<'a> {
    bam_chunk: BamChunk<'a>,
    header: HeaderView,
    filters: FiberFilters,
    cur_chunk: Vec<FiberseqData>,
}

impl<'a> FiberseqRecords<'a> {
    pub fn new(bam: &'a mut bam::Reader, filters: FiberFilters) -> Self {
        let header = bam.header().clone();
        let bam_recs = bam.records();
        let mut bam_chunk = BamChunk::new(bam_recs, None);
        bam_chunk.set_bit_flag_filter(filters.bit_flag);
        let cur_chunk: Vec<FiberseqData> = vec![];
        FiberseqRecords {
            bam_chunk,
            header,
            filters,
            cur_chunk,
        }
    }
}

impl Iterator for FiberseqRecords<'_> {
    type Item = FiberseqData;

    fn next(&mut self) -> Option<Self::Item> {
        // if we are out of data check for another chunk in the bam
        if self.cur_chunk.is_empty() {
            match self.bam_chunk.next() {
                Some(recs) => {
                    self.cur_chunk = FiberseqData::from_records(recs, &self.header, &self.filters);
                    // we will be popping from this list so we want to remove the first element first, not the last
                    self.cur_chunk.reverse();
                }
                None => return None,
            }
        }
        self.cur_chunk.pop()
    }
}

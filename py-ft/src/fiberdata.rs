use fibertools_rs::extract::FiberseqData;
use pyo3::iter::IterNextOutput;
use pyo3::prelude::*;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux, bam::Read};
//use rust_htslib::bam::Read;
//use std::io::Result;

#[pyclass]
pub struct PyFiberdata {
    // PB features
    #[pyo3(get, set)]
    pub ec: i64,
    // record features
    #[pyo3(get, set)]
    pub qname: String,
    #[pyo3(get, set)]
    pub qlen: i64,
    #[pyo3(get, set)]
    pub sam_flag: u16,
    #[pyo3(get, set)]
    pub rq: String,
    #[pyo3(get, set)]
    pub hp: String,
    #[pyo3(get, set)]
    pub seq: String,
    #[pyo3(get, set)]
    pub ct: String,
    #[pyo3(get, set)]
    pub start: i64,
    #[pyo3(get, set)]
    pub end: i64,
    #[pyo3(get, set)]
    pub strand: char,
    #[pyo3(get, set)]
    pub rg: String,
    // fiberseq features
    #[pyo3(get, set)]
    pub m6a: (Vec<i64>, Vec<i64>, Vec<u8>),
    #[pyo3(get, set)]
    pub msp: (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>),
    #[pyo3(get, set)]
    pub nuc: (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>),
}
#[pymethods]
impl PyFiberdata {
    pub fn __str__(&self) -> String {
        format!("{}\t{}\t{}\t{}", self.qname, self.ct, self.start, self.end,)
    }
}

fn new_py_fiberdata(fiber: &FiberseqData) -> PyFiberdata {
    // PB features
    let ec = fiber.ec.round() as i64;

    // record features
    let qname = std::str::from_utf8(fiber.record.qname()).unwrap();
    let qlen = fiber.record.seq_len() as i64;
    let sam_flag = fiber.record.flags();
    let rq = match fiber.get_rq() {
        Some(x) => format!("{}", x),
        None => ".".to_string(),
    };
    let hp = fiber.get_hp();
    let seq = String::from_utf8_lossy(&fiber.record.seq().as_bytes()).to_string();
    let (ct, start, end, strand): (&str, i64, i64, char) = if fiber.record.is_unmapped() {
        (".", 0, 0, '.')
    } else {
        let strand = if fiber.record.is_reverse() { '-' } else { '+' };
        (
            &fiber.target_name,
            fiber.record.reference_start(),
            fiber.record.reference_end(),
            strand,
        )
    };
    let rg = if let Ok(Aux::String(f)) = fiber.record.aux(b"RG") {
        log::trace!("{f}");
        f
    } else {
        "."
    };

    // fiberseq features
    let m6a = fiber.base_mods.m6a();
    let msp_starts = fiber.get_msp(false, true);
    let msp_lengths = fiber.get_msp(false, false);
    let ref_msp_starts = fiber.get_msp(true, true);
    let ref_msp_lengths = fiber.get_msp(true, false);
    let msp = (msp_starts, msp_lengths, ref_msp_starts, ref_msp_lengths);

    let nuc_starts = fiber.get_nuc(false, true);
    let nuc_lengths = fiber.get_nuc(false, false);
    let ref_nuc_starts = fiber.get_nuc(true, true);
    let ref_nuc_lengths = fiber.get_nuc(true, false);
    let nuc = (nuc_starts, nuc_lengths, ref_nuc_starts, ref_nuc_lengths);

    PyFiberdata {
        ec,
        qname: qname.to_string(),
        qlen,
        sam_flag,
        rq,
        hp,
        seq,
        ct: ct.to_string(),
        start,
        end,
        strand,
        rg: rg.to_string(),
        m6a,
        msp,
        nuc,
    }
}

#[pyclass]
pub struct FiberdataIter {
    count: usize,
    fiberdata: Vec<FiberseqData>,
}

#[pymethods]
impl FiberdataIter {
    #[new]
    pub fn new(f: &str, chrom: &str, start: i64, end: i64) -> Self {
        let mut bam = bam::IndexedReader::from_path(f).expect("unable to open indexed bam file");
        let header = bam::Header::from_template(bam.header());
        let head_view = bam::HeaderView::from_header(&header);
        bam.fetch((chrom, start, end))
            .expect("unable to fetch region");
        let records: Vec<bam::Record> = bam.records().map(|r| r.unwrap()).collect();
        let fiberdata = FiberseqData::from_records(&records, &head_view, 150);

        FiberdataIter {
            count: 0,
            fiberdata,
        }
    }

    fn __next__(&mut self) -> IterNextOutput<PyFiberdata, &'static str> {
        if self.count < self.fiberdata.len() {
            self.count += 1;
            // Given an instance `counter`, First five `next(counter)` calls yield 1, 2, 3, 4, 5.
            IterNextOutput::Yield(new_py_fiberdata(&self.fiberdata[self.count - 1]))
        } else {
            IterNextOutput::Return("Ended")
        }
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
}